"""
    R_p(ro, ri, kp)

Function that computes the pipe thermal resistance (radial conduction resistance around a
cylinder pipe). Valid for 1 pipe.
# Arguments
    - `ro`: Outer radius of the pipe [m]
    - `ri`: Inner radius of the pipe [m]
    - `kp`: Pipe thermal conductivity [W/mK]
# Output
    - `Rp`: Pipe conductive thermal resistance [mK/W]
# Reference
    - Bergman, T.L., Incropera, F.P.: Fundamentals of Heat and Mass Transfer, 7th edn. Wiley, New 
        York (2011)
    - Lamarche, L. (2023). Fundamentals of Geothermal Heat Pump Systems: Design and Application. 
        Springer Nature Switzerland.
"""
function resistance_pipe(ro::Real, ri::Real, kp::Real)
    return log(ro / ri) / (2 * π * kp)
end