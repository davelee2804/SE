#!/usr/bin/env bash

echo 'gdb --args ./kelvinhelmholtz 0 -helmholtz_ksp_type cg -helmholtz_pc_type jacobi -helmholtz_ksp_converged_reason -poisson_ksp_converged_reason -helmholtz_ksp_rtol 1.0e-16'
gdb --args ./kelvinhelmholtz 0 -helmholtz_ksp_type cg -helmholtz_pc_type jacobi -helmholtz_ksp_converged_reason -poisson_ksp_converged_reason -helmholtz_ksp_rtol 1.0e-16
