# NeuralNetwork/

In the reference toolbox this folder holds the neural-network representation of p(x) (`meshFlag = 2`).
It produces p̂(x) = Px + N(x)(x ⊙ x), one network output per costate component.

The soft pneumatic problem (ACNDC 2027) does **not** use a neural network. The paper fits one scalar,
positive-definite value function

    Vhat(y) = (y' P5 y) * exp(g(y)),   g(y) = sum_{1<=|alpha|<=3} c_alpha y^alpha   (55 coefficients)

so that the five costates are the gradient of a single function and Hess Vhat(0) = 2 P5 for any c.
That representation is implemented in `polybase/` (`f_value_eval.m`, `f_value_fit.m`) and selected with
`meshFlag = 0`. This folder is kept only so that the package has the same layout as `Stable-manifold-HJ-master`.
`P5_solution_check.m` stops with an error if `meshFlag` is not 0. The reference implementation is in
`Stable-manifold-HJ-master/NeuralNetwork/`.
