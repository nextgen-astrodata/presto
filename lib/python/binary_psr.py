import umath, Numeric, parfile, psr_utils
from psr_constants import *

def myasarray(a):
    if type(a) in [type(1.0),type(1L),type(1),type(1j)]:
        a = Numeric.asarray([a])
    if len(a) == 0:
        a = Numeric.asarray([a])
    return a

class binary_psr:
    """
    class binary_psr

        This class reads in a parfile (the only option for instantiation) of
            a binary pulsar.  It allows access the calculation of the mean,
            eccentric, and true anomalies, orbital position, radial velocity,
            and predicted spin period as a function of time.
            
    """

    def __init__(self, parfilenm):
        self.par = parfile.psr_par(parfilenm)
        if not self.par.__dict__.has_key('BINARY'):
            print "'%s' doesn't contain parameters for a binary pulsar!"
            return None
        else:
            if self.par.__dict__.has_key("TASC"):
                self.T0 = self.par.TASC
            else:
                self.T0 = self.par.T0
        self.PBsec = self.par.PB*SECPERDAY

    def calc_anoms(self, MJD):
        """
        calc_anoms(MJD):
            Return a tuple of the mean, eccentric, and true anomalies (all in radians)
                at the barycentric epoch MJD(s).
        """
        MJD = myasarray(MJD)
        difft = (MJD - self.T0)*SECPERDAY
        sec_since_peri = umath.fmod(difft, self.PBsec)
        if (sec_since_peri < 0.0): sec_since_peri += self.PBsec
        mean_anom = sec_since_peri/self.PBsec*TWOPI
        ecc_anom = self.eccentric_anomaly(mean_anom)
        true_anom = psr_utils.true_anomaly(ecc_anom, self.par.E)
        return (mean_anom, ecc_anom, true_anom)

    def most_recent_peri(self, MJD):
        """
        most_recent_peri(MJD):
            Return the MJD(s) of the most recent periastrons that occurred
            before the input MJD(s).
        """
        MJD = myasarray(MJD)
        difft = MJD - self.T0
        days_since_peri = umath.fmod(difft, self.par.PB)
        if (days_since_peri < 0.0): days_since_peri += self.par.PB
        return MJD - days_since_peri

    def eccentric_anomaly(self, mean_anomaly):
        """
        eccentric_anomaly(mean_anomaly):
            Return the eccentric anomaly in radians, given a set of mean_anomalies
            in radians.
        """
        ma = umath.fmod(mean_anomaly, TWOPI)
        ma = Numeric.where(ma < 0.0, ma+TWOPI, ma)
        eccentricity = self.par.E
        ecc_anom_old = ma
        ecc_anom = ma + eccentricity*umath.sin(ecc_anom_old)
        # This is a simple iteration to solve Kepler's Equation
        while (umath.maximum.reduce(umath.fabs(ecc_anom-ecc_anom_old)) > 1e-15):
            ecc_anom_old = ecc_anom[:]
            ecc_anom = ma + eccentricity*umath.sin(ecc_anom_old)
        return ecc_anom

    def calc_omega(self, MJD):
        """
        calc_omega(MJD):
            Return the argument of periastron (omega in radians) at time (or times) MJD(s).
        """
        MJD = myasarray(MJD)
        difft = (MJD - self.T0)*SECPERDAY
        if self.par.__dict__.has_key('OMDOT'):
            # Note:  This is an array
            return (self.par.OM + difft/SECPERJULYR*self.par.OMDOT)*DEGTORAD
        else:
            return self.par.OM*DEGTORAD

    def radial_velocity(self, MJD):
        """
        radial_velocity(MJD):
            Return the radial velocity of the pulsar (km/s) at the given MJD(s).
        """
        ma, ea, ta = self.calc_anoms(MJD)
        ws = self.calc_omega(MJD)
        c1 = TWOPI*self.par.A1/self.PBsec;
        c2 = umath.cos(ws)*umath.sqrt(1-self.par.E*self.par.E);
        sws = umath.sin(ws);
        cea = umath.cos(ea)
        return SOL/1000.0*c1*(c2*cea - sws*umath.sin(ea)) / (1.0 - self.par.E*cea)

    def doppler_period(self, MJD):
        """
        doppler_period(MJD):
            Return the observed pulse spin period in sec at the given MJD(s).
        """
        vs = self.radial_velocity(MJD)*1000.0 # m/s
        return self.par.P0*(1.0+vs/SOL)

    def position(self, MJD, inc=60.0):
        """
        position(MJD):
            Return the 'x' (along the LOS with + being towards us) and 'y' (in the
                plane of the sky with + being away from the line of nodes and -
                being in the direction of the line of nodes) positions of the
                pulsar with respect to the center of mass in units of lt-sec.
                (Note:  This places the observer at (+inf,0.0) and the line of nodes
                extending towards (0.0,-inf) with the pulsar orbiting (0.0,0.0)
                clockwise).  'inc' is the inclination of the orbit in degrees.
                MJD can be an array.  The return value is (xs, ys).
        """
        ma, ea, ta = self.calc_anoms(MJD)
        ws = self.calc_omega(MJD)
        orb_phs = ta + ws
        sini = umath.sin(inc*DEGTORAD)
        x = self.par.A1/sini
        r = x*(1.0-self.par.E*self.par.E)/(1.0+self.par.E*umath.cos(ta))
        return -r*umath.sin(orb_phs)*sini, -r*umath.cos(orb_phs)


if __name__=='__main__':
    from Pgplot import *
    
    psrA = binary_psr("0737A_Lyne_DD.par")
    times = psr_utils.span(0.0, psrA.par.PB, 1000) + psrA.par.T0
    rv = psrA.radial_velocity(times)
    plotxy(rv, (times-psrA.par.T0)*24,
           labx="Hours since Periastron", laby="Radial Velocity (km.s)")
    closeplot()
    print psrA.calc_anoms(52345.32476876)
    