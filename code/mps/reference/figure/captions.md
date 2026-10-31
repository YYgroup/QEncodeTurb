# Figure caption

Matrix product state (MPS) representations of geometrically generated initial velocity fields on $N^3$ grids, with $n=3\log_2N$ spatial qubits.
\textbf{a} Normalized Schmidt spectra of $u_x$ at the cuts requiring the largest single-cut rank for a relative $L^2$ error of $10^{-2}$. Lines and markers identify the grid size and cut, as indicated in the legend.
\textbf{b} Schmidt entropy $S_\ell$ at the relative scale cut $\ell/L$, where $L=\log_2N$. Colors and markers match \textbf{a}; the gray dotted line gives the entropy bound for $n=27$.
\textbf{c} Verified MPS bond dimensions normalized by the octal-encoding bound $\chi_{\max}=8^{\lfloor L/2\rfloor}$. Downward bars extend to the necessary ranks inferred from the Schmidt spectra and do not represent statistical uncertainty.
\textbf{d} Schmidt lower bound on the maximum componentwise relative $L^2$ error for the $512^3$ grid as a function of the bond cap $\chi$. Symbols give measured MPS reconstruction errors; colored horizontal lines mark the target errors.
\textbf{e,f} Gauge-reduced parameter counts and stored MPS scalars, respectively, normalized by the dense velocity representation with $3N^3$ real scalars. The horizontal reference in \textbf{c}, \textbf{e}, and \textbf{f} marks unity.
Each grid uses one realization. The population fractions of the first spinor component are $0.5$, $0.5$, $0.08$, $0.02$, and $0.006$ for increasing $N$.

# Table 1 caption

Resources of the velocity MPS at a componentwise relative $L^2$ error tolerance $\epsilon=10^{-2}$.
The spatial register contains $n$ qubits and represents an $N^3$ grid.
The octal scale encoding has $L=\log_2N$ sites and the rank bound $\chi_{\max}=8^{\lfloor L/2\rfloor}$.
The necessary bound $\chi_{\mathrm{lb}}$ follows from the Schmidt spectra, whereas $\chi_{\mathrm{built}}$ is attained by an explicitly constructed MPS with verified reconstruction error.
Both ranks are maxima over all cuts and velocity components.
$P_{\mathrm{NVPS}}$ excludes internal gauge freedom.
Parameter counts and storage sum the three velocity components.
Storage uses double precision, with $1\,\mathrm{MiB}=2^{20}$ bytes.
Each row represents one realization with the manuscript parameters.

# Table 2 caption

Resources of the generating spinor circuits.
The qubit count gives the spatial register size of each component circuit.
The parameter count includes shared amplitude angles, independent phase angles for both spinor components, and one population fraction.
Gate counts are summed over the two circuits, and the reported depth is their maximum.
Transpilation uses the $\{U,\mathrm{CX}\}$ basis at optimization level 3 without a hardware coupling constraint.
These counts cover spinor preparation before Gray decoding and Fourier transforms.
They exclude velocity evaluation, measurement, and time evolution.
