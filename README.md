# GroundHeatExchanger.jl

A Julia package for the thermal simulation of ground heat exchangers (GHEs). The package provides analytical ground models, temporal and spatial superposition techniques, and borehole thermal resistance calculations for single closed-loop borehole heat exchangers.

## Features

### Ground Models
Analytical solutions for thermal response of borehole heat exchangers:
- `ils` - Infinite line source (ILS) based on Ingersol (1948)
- `ics` - Infinite cylindrical source (ICS) based on Carlsaw & Jaeger (1959)
- `fls` - Finite line source (FLS) based on Claesson & Javed (2011)
- `mils` - Moving infinite line source (MILS) based on Pasquier & Lamarche (2022)
- `mfls` - Moving finite line source (MFLS) based on Guo et al. (2020)

All ground models support multiple dispatch for:
- Single or multiple time steps
- Single or multiple radii
- 2D arrays for borefield radius calculations
- Automatic type promotion and broadcasting

### Temporal Superposition
Fast convolution techniques for transient thermal response:
- `convolution` / `convolutionf` - FFT-based O(n log n) convolution for equally spaced time steps with stationary heat load profiles
- `convolution_ns` - Convolution for non-stationary heat load profiles with time-varying flowrates
- `state_vector` / `state_indices` - Helper functions for state transitions in non-stationary simulations
- `impulse_func` / `impulse_func_ns` - Impulse response functions for superposition

### Spatial Superposition (Borefield G-Functions)
Methods for computing thermal interactions in borehole fields:
- `bloc_matrix` - Block matrix formulation for g-function construction
- `successive_flux` - Iterative successive flux method for rapid g-function computation

### Borehole Thermal Resistance
Thermal resistance calculations for U-tube heat exchangers:
- `resistance_fluid` - Fluid-to-pipe convective resistance
- `resistance_pipe` - Pipe conductive resistance  
- `resistance_borehole_multipole` / `resistance_total_internal_multipole` - Multipole method for borehole thermal resistance (orders 0 and 1)
- `resistance_borehole_effective` - Effective (total) borehole thermal resistance
- `friction_factor_Colebrook_White` - Darcy friction factor calculation for pipe flow
- `Reynold`, `Prandtl`, `Nusselt` - Dimensionless flow characteristics

### Temperature Simulation
Fluid temperature response in borehole heat exchangers:
- `fluid_temperature` - Mean fluid temperature evolution over time
- `outlet_temperature` - Outlet fluid temperature
- `inlet_temperature` - Inlet fluid temperature

### Utility Functions
Helper functions for borefield design and analysis:
- `pchip_interpolation` - Piecewise Cubic Hermite Interpolating Polynomial for smooth data interpolation
- `set_nodes` - Logarithmic node positioning for transfer function computations
- `borefield_xy` - Generate rectangular borefield borehole coordinates
- `borefield_radius` - Calculate distances between boreholes in a borefield
- `water_ρ`, `water_cp`, `water_k`, `water_μ` - Water thermophysical properties as functions of temperature
- `head_loss_Darcy_Weisbach` - Pipe head loss calculations
- `heat_load_profile` - Sinusoidal heat load profile generation with seasonal modulation
- `GHE` - Pre-configured parameters for typical GHE systems

## Installation

The package is not yet registered. Install directly from the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/gabriel-dion/GroundHeatExchanger.jl")
```

Or in the Julia REPL package manager (`]`):

```
pkg> add https://github.com/gabriel-dion/GroundHeatExchanger.jl
```

## Dependencies

- [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl)
- [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl)
- [DSP.jl](https://github.com/JuliaDSP/DSP.jl)
- [PCHIPInterpolation.jl](https://github.com/gerlero/PCHIPInterpolation.jl)
- [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) (standard library)
- [CairoMakie](https://github.com/MakieOrg/Makie.jl)

## References

### Ground Thermal Response Models
- Ingersol, L. R. (1948). Theory of the ground pipe heat source for the heat pump. ASHVE Journal Section, Heating, Piping and Air Conditioning.
- Carslaw, H. S., & Jaeger, J. C. (1959). Conduction of heat in solids. Oxford: Clarendon Press, 1959, 2nd Ed.
- Claesson, J., & Javed, S. (2011). An analytical method to calculate borehole fluid temperatures for time-scales from minutes to decades. ASHRAE Transactions, 117(PART 2), 279–288.
- Guo, Y., Hu, X., Banks, J., & Liu, W. V. (2020). Considering buried depth in the moving finite line source model for vertical borehole heat exchangers—A new solution. Energy and Buildings, 214, 109859. https://doi.org/10.1016/j.enbuild.2020.109859
- Pasquier, P., & Lamarche, L. (2022). Analytic expressions for the moving infinite line source model. Geothermics, 103, 102413. https://doi.org/10.1016/j.geothermics.2022.102413

### Temporal and Spatial Superposition
- Marcotte, D., & Pasquier, P. (2008). Fast fluid and ground temperature computation for geothermal ground-loop heat exchanger systems. Geothermics, 37(6), 651–665. https://doi.org/10.1016/j.geothermics.2008.08.003
- Pasquier, P., & Marcotte, D. (2013). Efficient computation of heat flux signals to ensure the reproduction of prescribed temperatures at several interacting heat sources. Applied Thermal Engineering, 59(1–2), 515–526. https://doi.org/10.1016/j.applthermaleng.2013.06.018
- Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for efficient g-function construction. Renewable Energy, 121, 249–260. https://doi.org/10.1016/j.renene.2017.12.092
- Beaudry, G., Pasquier, P., & Nguyen, A. (2024). New formulations and experimental validation of non-stationary convolutions for the fast simulation of time-variant flowrates in ground heat exchangers. Science and Technology for the Built Environment, 30(3), 208–219. https://doi.org/10.1080/23744731.2023.2279468
- Nguyen, A., & Pasquier, P. (2021). A successive flux estimation method for rapid g-function construction of small to large-scale ground heat exchanger. Renewable Energy, 165, 359–368. https://doi.org/10.1016/j.renene.2020.10.074

### Borehole Thermal Resistance
- Claesson, J., & Hellström, G. (2011). Multipole method to calculate borehole thermal resistances in a borehole heat exchanger. Hvac&R Research, 17(6), 895–911.
- Javed, S., & Spitler, J. D. (2016). 3—Calculation of borehole thermal resistance. In S. J. Rees (Ed.), Advances in Ground-Source Heat Pump Systems (pp. 63–95). Woodhead Publishing. https://doi.org/10.1016/B978-0-08-100311-4.00003-0
- Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance calculation methods for grouted single U-tube ground heat exchangers. Applied Energy, 187, 790–806. https://doi.org/10.1016/j.apenergy.2016.11.079
- Claesson, J., & Javed, S. (2019). Explicit multipole formulas and thermal network models for calculating thermal resistances of double U-pipe borehole heat exchangers. Science and Technology for the Built Environment, 25(8), 980–992. https://doi.org/10.1080/23744731.2019.1620565

### Foundational References
- Hellström, Göran. 1991. "Ground Heat Storage : Thermal Analyses of Duct Storage Systems." http://www.lunduniversity.lu.se/o.o.i.s?id=24732&postid=2536279.
- Lamarche, L. (2023). Fundamentals of Geothermal Heat Pump Systems: Design and Application. Springer Nature Switzerland.
- Carslaw, H. S., & Jaeger, J. C. (1959). Conduction of heat in solids. Oxford: Clarendon Press, 1959, 2nd Ed.
- Claesson, J., & Javed, S. (2011). An analytical method to calculate borehole fluid temperatures for time-scales from minutes to decades. ASHRAE Transactions, 117(PART 2), 279–288.
- Pasquier, P., & Lamarche, L. (2022). Analytic expressions for the moving infinite line source model. Geothermics, 103, 102413. https://doi.org/10.1016/j.geothermics.2022.102413
- Guo, Y., Hu, X., Banks, J., & Liu, W. V. (2020). Considering buried depth in the moving finite line source model for vertical borehole heat exchangers—A new solution. Energy and Buildings, 214, 109859. https://doi.org/10.1016/j.enbuild.2020.109859
- Javed, S., & Spitler, J. (2017). Accuracy of borehole thermal resistance calculation methods for grouted single U-tube ground heat exchangers. Applied Energy, 187, 790–806. https://doi.org/10.1016/j.apenergy.2016.11.079
- Hellström, Göran. 1991. “Ground Heat Storage : Thermal Analyses of Duct Storage Systems.” http://www.lunduniversity.lu.se/o.o.i.s?id=24732&postid=2536279.
- Lamarche, L. (2023). Fundamentals of Geothermal Heat Pump Systems: Design and Application. Springer Nature Switzerland.
- Claesson, J., & Hellström, G. (2011). Multipole method to calculate borehole thermal resistances in a borehole heat exchanger. Hvac&R Research, 17(6), 895–911.
- Javed, S., & Spitler, J. D. (2016). 3—Calculation of borehole thermal resistance. In S. J. Rees (Ed.), Advances in Ground-Source Heat Pump Systems (pp. 63–95). Woodhead Publishing. https://doi.org/10.1016/B978-0-08-100311-4.00003-0
- Marcotte, D., & Pasquier, P. (2008). Fast fluid and ground temperature computation for geothermal ground-loop heat exchanger systems. Geothermics, 37(6), 651–665. https://doi.org/10.1016/j.geothermics.2008.08.003
- Pasquier, P., & Marcotte, D. (2013). Efficient computation of heat flux signals to ensure the reproduction of prescribed temperatures at several interacting heat sources. Applied Thermal Engineering, 59(1–2), 515–526. https://doi.org/10.1016/j.applthermaleng.2013.06.018
- Beaudry, G., Pasquier, P., & Nguyen, A. (2024). New formulations and experimental validation of non-stationary convolutions for the fast simulation of time-variant flowrates in ground heat exchangers. Science and Technology for the Built Environment, 30(3), 208–219. https://doi.org/10.1080/23744731.2023.2279468
- Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for efficient g-function construction. Renewable Energy, 121, 249–260. https://doi.org/10.1016/j.renene.2017.12.092
- Nguyen, A., & Pasquier, P. (2021). A successive flux estimation method for rapid g-function construction of small to large-scale ground heat exchanger. Renewable Energy, 165, 359–368. https://doi.org/10.1016/j.renene.2020.10.074
- Claesson, J., & Javed, S. (2019). Explicit multipole formulas and thermal network models for calculating thermal resistances of double U-pipe borehole heat exchangers. Science and Technology for the Built Environment, 25(8), 980–992. https://doi.org/10.1080/23744731.2019.1620565

## License

See [LICENSE.txt](LICENSE.txt).
