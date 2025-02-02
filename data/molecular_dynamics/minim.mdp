; md.mdp - Production MD run parameters
integrator      = md         ; Molecular Dynamics integrator
nsteps          = 50000      ; Number of MD steps (adjust as necessary)
dt              = 0.001      ; Time step in ps
nstxout         = 1000       ; Save coordinates every 1000 steps
nstvout         = 1000       ; Save velocities every 1000 steps
nstenergy       = 1000       ; Save energies every 1000 steps
nstlog          = 1000       ; Write to log file every 1000 steps
continuation    = yes
constraints     = all-bonds
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = System
tau_t           = 0.1
ref_t           = 300

; Pressure coupling
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5