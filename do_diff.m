syms u
syms v
f = .5*u^3*(1-v)^3 / ((1-v)^3*(1-u)^3+.5*u^3*(1-v)^3)
fu   = simplify(diff(f  , u))
fv   = simplify(diff(f  , v))

fuu  = simplify(diff(fu , u))
fuv  = simplify(diff(fu , v))
fvv  = simplify(diff(fv , v))

fuuu = simplify(diff(fuu, u))
fuuv = simplify(diff(fuu, v))
fuvv = simplify(diff(fuv, v))
fvvv = simplify(diff(fvv, v))
