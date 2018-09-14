# ragnar_imsrg
In-Medium Similarity Renormalization Group software for nuclear structure calculations. Written in C++, with Python bindings.
It is capable of performing Hartree-Fock, single-reference IM-SRG, and valence-space IM-SRG calculations.
If you don't know what these things are, you should check out
http://link.aps.org/doi/10.1103/PhysRevLett.106.222502
http://www.sciencedirect.com/science/article/pii/S0370157315005414

The IMSRG code in the src directory is copyrighted under the GNU Public License (http://www.gnu.org/licenses/gpl.html).
This code also uses the Armadillo library, which is covered under the Mozilla Public License (https://www.mozilla.org/en-US/MPL/2.0/), and the BOOST library, which is covered under the BOOST Software License (http://www.boost.org/users/license.html)


Credits:
	* Most of this code was written by Ragnar Stroberg, who benefitted greatly from the codes of Koshiroh Tsukiyama, Heiko Hergert, and Nathan Parzuchowski.
	* Implementation of operators for neutrinoless double beta decay was done by Charlie Payne.
	* Some tweaks to the 3rd order MBPT code were contributed by Johannes Simonis.
