from sympy import *
import numpy as np
import matplotlib.pyplot as plt

"""
Trouver la contrainte dans les cas suivant:

1. Configuration en vol de croisière (équilibre).
2. Configuration en équilibre avec un facteur de charge de 3g (poids x3).
3. Avion au sol soutenu par son train d’atterrissage.
"""


class AirPlane:
    def __init__(self):
        self.chargeFactor = 1
        self.onGround = False
        self.resolution = 1000

        self.wingLength = 150
        self.totalWeight = 0
        self.nervures = np.array([[25, 50, 80.8, 111.5, 150], [7.13, 5.89, 6.86, 3.35, 3.21]])
        self.wheel = [25, 110]
        self.bl = 1.648
        self.tl = 0.16

        self.pAlum = 0.101
        self.pFuel = 0.0303
        self.beamHeightFunc = None
        self.beamAreaFunc = None
        self.shear = []
        self.moment = [0]

    def getWingShearAndMoment(self):
        if self.onGround:
            self.nervures[1][0] -= 1815.96
        self.getTotalWeight()

        A = np.linspace(150, 0, self.resolution)[1:]

        """Calcul de l'effort tranchant"""

        for i, a in enumerate(A):
            if self.onGround:
                shear = self.getWingSectionWeight(a)
            else:
                shear = self.getWingSectionWeight(a) - self.getAeroLoad(a)
                shear *= self.chargeFactor

            self.shear.append(-1*shear)

        fig, (ax1, ax2) = plt.subplots(2, figsize=(7, 6))
        ax1.plot(A, self.shear)
        ax1.set_ylabel("Effort Tranchant [lb]", fontsize=13)

        """Calcul du moment"""

        X = A
        # X = list(reversed(A))
        # self.shear = list(reversed(self.shear))
        Xstep = 150 / (self.resolution-1)

        for y in self.shear:
            self.moment.append(y*Xstep + self.moment[-1])

        ax2.plot(X, self.moment[1:])
        ax2.set_ylabel("Moment [lb$\\times$in]", fontsize=13)

        ax2.set_xlabel("Distance sur l'aile [pouce]", fontsize=13)

        plt.show()

    def getBeamStrain(self):
        X = np.linspace(150, 0, self.resolution)[1:]
        x = symbols("x")

        inertiaFunc = (0.0018394 * 0.125 * (((-20.6/150) * x) + 41.2)**3) + 2*(((1/12) * self.bl * (self.beamHeightFunc**3)) - ((1/12) * (self.bl - self.tl) * (self.beamHeightFunc - 2*self.tl)**3))
        QFunc = ((self.bl * self.tl * (self.beamHeightFunc/2 + self.tl/2)) + (self.tl * (self.beamHeightFunc/2 - self.tl) * ((self.beamHeightFunc/2 - self.tl/2)/2))) * 2

        Q = []
        height = []
        inertia = []
        for d in X:
            Q.append(QFunc.evalf(subs={x: d}))
            height.append(self.beamHeightFunc.evalf(subs={x: d}))
            inertia.append(inertiaFunc.evalf(subs={x: d}))

        normalStrain = ((np.array(self.moment[1:]) * (np.array(height)/2)) / (np.array(inertia))) / 1000
        shearStrain = ((np.array(self.shear) * np.array(Q)) / (np.array(inertia) * self.tl)) / 1000

        print("Max Normal Strain = {} || SF = {}".format(round(np.max(normalStrain), 3), round(60/np.max(normalStrain), 3)))
        print("Max Shear  Strain = {} || SF = {}".format(round(np.max(shearStrain), 3), round(25/np.max(shearStrain), 3)))

        fig, (ax1, ax2) = plt.subplots(2, figsize=(7, 6))
        ax1.plot(X, shearStrain)
        ax2.plot(X, normalStrain)
        ax2.set_ylabel("Contrainte\nnormale [ksi]", fontsize=13)
        ax1.set_ylabel("Contrainte\nde cisaillement [ksi]", fontsize=13)
        ax2.set_xlabel("Distance sur l'aile [pouce]", fontsize=13)
        plt.show()

    def getTotalWeight(self):
        mainWeight = 5200
        wingWeight = self.getWingSectionWeight()
        self.totalWeight = 2 * wingWeight + mainWeight
        print("wing/Total = {} / {}".format(wingWeight, self.totalWeight))

    def getWingSectionWeight(self, a=0, b=150):
        fuelWeight = self.getFuelVolume(a, b) * self.pFuel
        beamWeight = 2 * self.getBeamVolume(a, b) * self.pAlum
        coatWeight = self.getCoatingVolume(a, b) * self.pAlum
        nervWeight = sum(self.nervures[1][np.where((self.nervures[0] >= a) & (self.nervures[0] <= b))[0]])
        wheelWeight = self.wheel[1] if a <= self.wheel[0] else 0

        weight = sum([fuelWeight, beamWeight, coatWeight, nervWeight, wheelWeight])

        return weight

    def getFuelVolume(self, a=0, b=150):
        x = symbols('x')
        beamDistance = (((7.2-14.4)/150) * x) + 14.4
        beamHeight = (((1.648-3.256)/150) * x) + 3.256

        beamAreaFunc = ((beamHeight-2*self.tl)*self.tl)+2*self.bl*self.tl

        areaFunc = beamDistance * beamHeight - 2*beamAreaFunc

        fuelVolume = integrate(areaFunc, (x, a, b))
        # meanX = integrate(x*areaFunc, (x, a, b)) / fuelVolume

        return fuelVolume  # , meanX

    def getBeamVolume(self, a=0, b=150):
        x = symbols('x')
        self.beamHeightFunc = (((1.648-3.256)/150) * x) + 3.256
        self.beamAreaFunc = ((self.beamHeightFunc-2*self.tl)*self.tl)+2*self.bl*self.tl

        beamVolume = integrate(self.beamAreaFunc, (x, a, b))
        # meanX = integrate(x*beamAreaFunc, (x, a, b)) / beamVolume

        return beamVolume  # , meanX

    def getAeroLoad(self, a=0, b=150):
        x = symbols('x')
        liftFunc = ((9 * self.totalWeight) / (16 * self.wingLength)) * (1 - (x / self.wingLength)**8)# fonction d prof
        aeroLoad = integrate(liftFunc, (x, a, b))
        # meanX = integrate(x*liftFunc, (x, a, b)) / aeroLoad

        return aeroLoad  # , meanX

    @staticmethod
    def getCoatingVolume(a=0, b=150):
        x = symbols('x')
        areaFunc = 2.024 * 0.125 * ((-20.6*x/150) + 41.2)

        coatVolume = integrate(areaFunc, (x, a, b))
        # meanX = integrate(x*areaFunc, (x, a, b)) / coatVolume

        return coatVolume  # , meanX


plane = AirPlane()

plane.getWingShearAndMoment()
plane.getBeamStrain()
