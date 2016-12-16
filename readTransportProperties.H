
    // Heat capacity
    dimensionedScalar Cp0("Cp0", dimArea/dimTime/dimTime/dimTemperature, laminarTransport);

    // Fluid density
    dimensionedScalar rho0("rho0", dimMass/dimVolume, laminarTransport);


    // Thermal expansion coefficient [1/K]
    dimensionedScalar beta
    (
        "beta",
        dimless/dimTemperature,
        laminarTransport
    );

    // Reference temperature [K]
    dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

    // Laminar Prandtl number
    dimensionedScalar Pr("Pr", dimless, laminarTransport);

    // Turbulent Prandtl number
    dimensionedScalar Prt("Prt", dimless, laminarTransport);
