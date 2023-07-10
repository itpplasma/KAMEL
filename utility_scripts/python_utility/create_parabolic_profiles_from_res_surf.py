import numpy as np

def create_parabolic_profiles_from_res_surf(path, q0, n0, Te0, Ti0, Vz0, Er0, Vth0, m_mode, n_mode, rmin, rmax, num, a, const=''):
	""" Create parabolic profiles for fixed density and electron
		temperature values at the rational surface. """
	r = np.linspace(rmin, rmax, num)
	q = -(1.05 + q0 * (r/a)**2)
	rres = np.interp(-m_mode/n_mode, r, q)
	
	fac_par = 1 - (r/a)**2
	n   = n0   * fac_par / (1-(rres/a)**2)
	Te  = Te0  * fac_par / (1-(rres/a)**2)
	Ti  = Ti0  * fac_par
	Vz  = Vz0  * fac_par
	Er  = Er0  * fac_par
	Vth = Vth0 * fac_par

	#if not const == '':
	#	print(f'profile {const} is set constant')
	np.savetxt(path + 'q.dat', np.array((r,q)).transpose())
	np.savetxt(path + 'Te.dat', np.array((r, Te)).transpose())
	np.savetxt(path + 'Ti.dat', np.array((r, Ti)).transpose())
	np.savetxt(path + 'n.dat', np.array((r, n)).transpose())
	np.savetxt(path + 'Vz.dat', np.array((r, Vz)).transpose())
	np.savetxt(path + 'Er.dat', np.array((r, Er)).transpose())
	np.savetxt(path + 'Vth.dat', np.array((r, Vth)).transpose())

