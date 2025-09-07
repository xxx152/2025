Consider：
$$
\omega_{n+1}(x) = \prod_{i=0}^n (x - x_i)
$$

$$
\omega_{n+1}(x) = (x - x_i) \prod_{j \neq i} (x - x_j)
$$

Taking the derivative with respect to $x$：
$$
\omega_{n+1}'(x) 
= \prod_{j \neq i} (x - x_j) + (x - x_i)\cdot \frac{d}{dx}\Bigg(\prod_{j \neq i} (x - x_j)\Bigg)
$$

At $x = x_i$, the second term  $(x - x_i)=0$. Thus:
$$
\omega_{n+1}'(x_i) = \prod_{j \neq i} (x_i - x_j).
$$

Therefore：
$$
w_i = \prod_{j \neq i} \frac{1}{x_i - x_j} 
= \frac{1}{\omega_{n+1}'(x_i)}.
$$
