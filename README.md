# GK
Space used to further understand Generalised Kelvin form of the creep compliance.

Given a parameterisation of the relaxation spectrum $X(\tau)$, it is possible to construct the time-domain creep compliance via the expression

$$ J(t) = J_U + \int X(\tau) (1 - e^{-t/\tau}) d (\text{ln} \tau) $$

Analogously to placing $N$ Maxwell (spring and dashpot in series) elements in parallel to create a Generalised Maxwell Model, we can place $N$ Kelvin (spring and dashpot in parallel) elements in series to create a Generalised Kelvin (GK) model. The corresponding mathematical form for the creep compliance is

$$ J(t) = J_U + \Sigma_i^N J_i (1 - e^{-t/\tau_i}) $$

Knowing that in the limit $N \rightarrow \infty$, the GK model approaches the same behaviour as a medium representable by a continuum relaxation spectrum, we can relate these two expressions by thinking of the second as an approximation of the first with finite $N$. In this case, we see that each component

$$ J_i = X(\tau_i) d (\text{ln} \tau_i) $$

where the interval $d (\text{ln} \tau_i)$ is the interval between successive $\tau_i$ in natural logarithm space. Since $X(\tau_i)$ is directly calculable from experimental parameterisations of anelasticity such as Yamauchi and Takei (2016), we can directly find the homologous-temperature- and Maxwell timescale- dependent coefficients $J_i$. The exact same approach can be applied to find $M(t)$ if the spectrum $Y(\tau)$ is known. In this situation, we can find the components according to

$$ M_i = Y(\tau_i) d (\text{ln} \tau_i) $$

Such values of $M_i$ correspond to the coefficients in the Generalised Maxwell approximation given by

$$ M(t) = M_U - \Sigma_i^N M_i (1 - e^{-t/\tau_i}) $$
