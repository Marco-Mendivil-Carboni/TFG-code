# Elastic potential

$$
V=\sum_{i=1}^{N-1}\frac{1}{2}k_e(\sqrt{(x_{i+1}-x_i)^2+(y_{i+1}-y_i)^2+(z_{i+1}-z_i)^2}-b)^2
$$

$$
\dot{p_{x_1}}=F_{x_1}=-\frac{\partial V}{\partial x_1}=k_1(x_2-x_1)
$$

$$
\dot{p_{x_i}}=F_{x_i}=-\frac{\partial V}{\partial x_i}=-k_{i-1}(x_i-x_{i-1})+k_i(x_{i+1}-x_i)
$$

$$
\dot{p_{x_N}}=F_{x_N}=-\frac{\partial V}{\partial x_N}=-k_{N-1}(x_N-x_{N-1})
$$

$$
k_{i}=\frac{k_e(\sqrt{(x_{i+1}-x_i)^2+(y_{i+1}-y_i)^2+(z_{i+1}-z_i)^2}-b)}{\sqrt{(x_{i+1}-x_i)^2+(y_{i+1}-y_i)^2+(z_{i+1}-z_i)^2}}
$$

# Bending potential

$$
V=\sum_{i=1}^{N-1} k_b (1-\cos(\theta_i))
$$

Let's define

$$
\bar{b}_i = \bar{r}_{i+1}-\bar{r}_i \quad , \quad d_i = \left\|\bar{b}_i\right\|
$$

then

$$
\cos(\theta_{i-1}) = \frac{\bar{b}_{i-1}\cdot\bar{b}_{i-2}}{d_{i-1}d_{i-2}} \quad , \quad \cos(\theta_i) = \frac{\bar{b}_i\cdot\bar{b}_{i-1}}{d_id_{i-1}} \quad , \quad \cos(\theta_{i+1}) = \frac{\bar{b}_{i+1}\cdot\bar{b}_{i}}{d_{i+1}d_{i}}
$$

the force will be

$$
-\frac{\partial V}{\partial r_i^j} = k_b \frac{\partial}{\partial r_i^j}\left(\cos(\theta_{i-1})+\cos(\theta_i)+\cos(\theta_{i+1})\right)
$$
$$
\frac{\partial}{\partial r_i^j}\cos(\theta_{i-1}) = \frac{1}{d_{i-1}d_{i-2}}\left(b_{i-2}^j-(\bar{b}_{i-1}\cdot\bar{b}_{i-2})b_{i-1}^j/d_{i-1}^2\right)  
$$
$$
\frac{\partial}{\partial r_i^j}\cos(\theta_i) = \frac{1}{d_id_{i-1}}\left(b_i^j-b_{i-1}^j+(\bar{b}_i\cdot\bar{b}_{i-1})\left(b_i^j/d_i^2-b_{i-1}^j/d_{i-1}^2\right)\right)  
$$
$$
\frac{\partial}{\partial r_i^j}\cos(\theta_{i+1}) = \frac{1}{d_{i+1}d_{i}}\left(-b_{i+1}^j+(\bar{b}_{i+1}\cdot\bar{b}_{i})b_{i}^j/d_{i}^2\right)  
$$

# Lennard-Jones potential (excluded volume)

$$
V=4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^{6}\right]
$$

where

$$
r=\sqrt{\sum_{j=0}^{2}(r_i^j-r_n^j)^{2}}
$$

thus

$$
-\frac{\partial V}{\partial r_i^j} = 4\varepsilon\left[ 12\frac{\sigma^{12}}{r^{14}}-6\frac{\sigma^{6}}{r^{8}}\right] (r_i^j-r_n^j)
$$
