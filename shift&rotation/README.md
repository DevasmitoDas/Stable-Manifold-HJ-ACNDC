# shift&rotation/

In the reference toolbox this folder holds the variational "shift / rotate and reform" tools.
They pick a computed trajectory, solve the variational BVP along it and re-converge a shifted or rotated
copy to fill gaps in the x-space coverage.

The ACNDC 2027 study does not use them. Its manifold anchors are the 18 accepted trajectories of P4
(radius 2e-4, no negative-time extension). Away from the origin, the value function is shaped by HJB
collocation and the guarded selector (paper Sec. 4). The paper names, as future work, extending the anchors along the
characteristics away from the origin. That would be the role of this folder, which holds only this README so that the package has the same layout as
`Stable-manifold-HJ-master`. The reference implementation is in `Stable-manifold-HJ-master/shift&rotation/`.
