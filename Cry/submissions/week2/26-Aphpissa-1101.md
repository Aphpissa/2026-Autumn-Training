# dp泄露攻击

*最燃尽的一集*

## 一.常规dp攻击

### 1.什么是dp？

dp是**CRT-RSA**中的量，是私钥的一部分

### 2.什么是CRT-RSA？

**CRT-RSA**是根据中国剩余定理（*CRT*）优化的RSA算法，这个算法可以让RSA的解密过程提高3~4倍

具体工作原理如下

我们设：

$$
dp = d \\bmod (p-1)
$$

$$
dq = d \\bmod (q-1)
$$

$$
qinv = p^{-1} \\bmod q
$$

这个环节就相当于把d拆成了两部分

原理就在于根据费马小定理

$$
c^{p=1}\equiv 1\pmod p
$$

而

$$
d\equiv dp\pmod p
$$

$$
q\equiv dq\pmod q
$$

所以

$$
c^{dp}\equiv c^d\pmod p
$$

同理，有

$$
c^{dq}\equiv c^d\pmod q
$$

所以

$$
m_p\equiv c^d\equiv m\pmod p
$$

$$
m_q\equiv c^d\equiv m\pmod q
$$

因而可以用中国剩余定理求解

### 3.持有私钥的解密过程

第一步，先算影子，即先对 $c$ 取 $dp,dq$ 的模

$$
m_p=c^{dp}\pmod p
$$

$$
m_q=c^{dq}\pmod q
$$

那怎么求解这个m呢

我们肯定会想到标准解法

<img width="779" height="344" alt="image" src="https://github.com/user-attachments/assets/a0174bda-7448-429d-b972-d05d2d43da16" />

这个解法肯定可行，但是计算量偏大，要存储的量偏多

这个时候我们就需要用到中国剩余定理的Garner解

就是先配一个满足第一个同余方程的解，再去配另一个同余方程的解

我们可以配出

$$
m=m_p+hp
$$

显然，对于这个 $m$ 来说，它永远满足第一个同余方程

那么接下来我们就要配出基于这个解的第二个同余方程的解

即找出 $h$ 使得

$$
m_q\equiv m_p+hp\pmod q
$$

显然可以找出来

直接进行移项

$$
m_q-m_p\equiv hp\pmod q
$$

然后乘上 $qinv$ 可得

$$
qinv(m_q-m_p)\equiv h\pmod q
$$

这就把 $h$ 给求出来了

至此，解密完成

$$
m=m_p+hp
$$

其中

$$
h\equiv qinv(m_q-m_p)\pmod q
$$

由于把 $d$ 拆成两份的缘故，数字变小，所以这样可以更快解密，比常规快出3~4倍

### 4.dp泄露

在CRT-RSA中， $dp$ （或 $dq$ ）都是需要保护的对象，不然如果 $dp$ （或 $dq$ ）泄露了，那么就和 $d$ 泄露没有什么区别

因为由

$$
ed\equiv 1\bmod (p-1)
$$

根据 $dp$ 的定义可知、

$$
dp=d\bmod (p-1)
$$

相乘可得

$$
eddp\equiv d\bmod (p-1)
$$

消去 $d$ 可得

$$
edp-1\equiv 0\bmod(p-1)
$$

也就是说， $edp-1$ 是 $p-1$ 的 $k$ 倍

也就可以用for循环在 $(0,e)$ 的范围内爆破出来！

### 5.实操环节

送上题目

```python

import gmpy2
import libnum
import uuid

flag = "flag{" + str(uuid.uuid4()) + "}"
print(flag)
m = libnum.s2n(flag)
p = libnum.generate_prime(1024)
q = libnum.generate_prime(1024)
e = 65537
n = p * q
phi = (p - 1) * (q - 1)
d = gmpy2.invert(e, phi)
print("d=",d)
dp = d % (p - 1)
c = pow(m, e, n)
print("n=", n)
print("e=", e)
print("c=", c)
print("dp=", dp)

```

由于这个 $k$ 解出来的数有可能不是 $n$ 的因子，因此需要进行判断

我的解法也贴在这里

```python

from Crypto.Util.number import long_to_bytes
import gmpy2

#依旧生成数

for k in range(1,e):
    if (e * dp - 1)%k==0:
        p=((e*dp-1)//k)+1
        if n%p!=0:
            continue
        q=n//p
        phi=(p-1)*(q-1)
        d=gmpy2.invert(e,phi)
        print(d)
        m=pow(c,d,n)
        print(long_to_bytes(int(m)))
        break
else:
    print("No solution")
```

运行结果符合格式

<img width="328" height="32" alt="image" src="https://github.com/user-attachments/assets/eb288505-b00a-4b55-ae47-2a81f29d3705" />

代码和答案差不多

就是要爆破一个倍数，且满足解出来的数是 $n$ 的因子

## 二.e过大

### 1.dp泄露的另解

当 $e$ 过大的时候，直接遍历在 $(0,e)$ 内的所有数肯定不现实

属于理论可行但CPU不接受理论的那种

因此，我们也可以用费马小定理来找一个其他解法

找一个数 $a$ ，且 $a$ 与 $p$ 互素

这时候就可以根据费马小定理得到

$$
a^{edp}-a\equiv 0\pmod p
$$

放在代码里，就是生成一个质数，代入上面那个式子以后就去和 $n$ 求公因数即可解出 $p$

但是有一个很复杂的点

电脑不可能运行一串超大串的幂运算， $a^{edp}$ 就直接算不出来了

解决方法是加入快速幂pow

快速幂pow可以一直把幂控制在小于10的范围内

且根据欧几里得算法，求公约数的时候如果其中一个数减去或加上另一个数的倍数，那么公约数不变

就可以解决诸如此类的问题

### 2.实操环节

题目

```python

import gmpy2
import libnum
import uuid

flag = "flag{" + str(uuid.uuid4()) + "}"
m = libnum.s2n(flag)
p = libnum.generate_prime(1024)
q = libnum.generate_prime(1024)
e = libnum.generate_prime(128)
n = p * q
phi = (p - 1) * (q - 1)
d = gmpy2.invert(e, phi)

dp = d % (p - 1)
c = pow(m, e, n)
print("n=", n)
print("e=", e)
print("c=", c)
print("dp=", dp)

```

这里采用libnum工具生成 $a$ 来保证 $p$ 与 $a$ 互素

```python

from Crypto.Util.number import getPrime,long_to_bytes
import gmpy2
import libnum

#生成的数

a = libnum.generate_prime(256)
print(a)
mid=pow(a,e*dp,n)-a
p=gmpy2.gcd(n,mid)
q=n//p
phi=(p-1)*(q-1)
d=gmpy2.invert(e,phi)
m=pow(c,d,n)
print(long_to_bytes(int(m)))

```

结果也是成功运行

<img width="286" height="34" alt="image" src="https://github.com/user-attachments/assets/5f088d76-bf9d-4940-861d-31aff2f6e84d" />

就是按照上述思路来的代码

只要有思路，剩下的都是翻译问题

## 三.dp高位泄露

### 1.什么是高位泄露

高位泄露，就是 $dp$ 的前几位泄露了，但是没有泄露完全

与之相对的，有低位泄露，也就是后几位泄露了

### 2.高位泄露怎么攻击

如果为高位泄露，那么我们就需要找到低位来填充整个 $dp$ 

根据上面的求解方法，
