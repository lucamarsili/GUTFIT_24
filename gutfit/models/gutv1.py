from gutfit.model import Model
import cmath

from numba import jit

@jit
def t11(x_):
    return x_[0]*(
            complex(-0.00014961882506530676,7.299895968258247e-7)/( x_[1]*(-1 + complex(0,1)*x_[3] + x_[2]) ) +
            complex(0.00002336*(-3 + complex(0,1)*x_[3] + x_[2]))/(-1 + complex(0,1)*x_[3] + x_[2]) +
            (x_[1]*(
                 (complex(2.4079600804265915e-9,1.1494748285782534e-17) - complex(2.409610269404143e-7,-8.514396563360973e-12)*x_[1] +
                     2.92e-6*x_[1]**2)/(
                         complex(-0.0014272717944412566,-1.793739297851695e-21) + complex(0.16855723761492827,0.005951877656781462)*x_[1] -
            complex(3.2849786279185222,-0.015627693218553906)*x_[1]**2 +
            1.*x_[1]**3))*(-1 + complex(0,1)*x_[3] +
                x_[2]))/(complex(0,1)*x_[3] + x_[2])
            )

@jit
def t12(x_):
    return x_[0]*(
            complex(-0.00032184486273232867,0.000011734522937336278)/(x_[1]*(-1 + complex(0,1)*x_[3] + x_[2])) +
            (x_[1]*(
                (complex(2.5823323634850127e-7,-1.2436436176592075e-10) - complex(0.000020115303920770545,-7.334076835835173e-7)*x_[1])/(
                    complex(0.0014272717944412566,1.793739297851695e-21) - complex(0.16855723761492827,0.005951877656781462)*x_[1] +
                complex(3.2849786279185222,-0.015627693218553906)*x_[1]**2 -
                1.*x_[1]**3))*(-1 + complex(0,1)*x_[3] + x_[2]))/(complex(0,1)*x_[3] + x_[2]))

@jit
def t13(x_):
    return x_[0]*(
            complex(-0.00011250381136069745,0.00032719113650081136)/(x_[1]*(-1 + complex(0,1)*x_[3] + x_[2])) +
            (x_[1]*(
                (complex(2.9669022224442643e-6,1.2949030596424693e-6) + complex(7.031488210043591e-6,-0.00002044944603130071)*x_[1])/(
                    complex(-0.0014272717944412566,-1.793739297851695e-21) +  complex(0.16855723761492827,0.005951877656781462)*x_[1] -
                    complex(3.2849786279185222,-0.015627693218553906)*x_[1]**2 +
                    1.*x_[1]**3))*(-1 + complex(0,1)*x_[3] + x_[2])
                )/(complex(0,1)*x_[3] + x_[2])
                )

@jit
def t22(x_):
    return x_[0]*(
            complex(-0.0015906758349636073,6.689733493850993e-8)/(x_[1]*(-1 + complex(0,1)*x_[3] + x_[2])) +
            (0.01144*(-3 + complex(0,1)*x_[3] + x_[2]))/(-1 + complex(0,1)*x_[3] +
        x_[2]) + (x_[1]*( (complex(0.0000598643192345929,-2.665466091266395e-8) -
            complex(0.004598102198238262,-0.000022343420219098426)*x_[1] +
            0.0014299999999999998*x_[1]**2)/(complex(-0.0014272717944412566,-1.793739297851695e-21)
                + complex(0.16855723761492827,0.005951877656781462)*x_[1] -
                complex(3.2849786279185222,-0.015627693218553906)*x_[1]**2 +
                1.*x_[1]**3))*(-1 + complex(0,1)*x_[3] + x_[2]))/(complex(0,1)*x_[3] + x_[2]))

@jit
def t23(x_):
    return x_[0]*(
            complex(-0.003930681666107606,9.314245518884659e-7)/(x_[1]*(-1 + complex(0,1)*x_[3] + x_[2])) +
            (x_[1]*((complex(0.0007434368991203491,0.00013861334759940198) - complex(0.0002456676041317253,-5.821403449302912e-8)*x_[1])/(
                complex(0.0014272717944412566,1.793739297851695e-21) - complex(0.16855723761492827,0.005951877656781462)*x_[1] +
                complex(3.2849786279185222,-0.015627693218553906)*x_[1]**2 - 1.*x_[1]**3))*(-1 + complex(0,1)*x_[3] + x_[2]))/(complex(0,1)*x_[3] + x_[2]))

@jit
def t33(x_):
    return x_[0]*(
            complex(-0.1110579400167097,-6.788494192461779e-8)/(x_[1]*(-1 + complex(0,1)*x_[3] + x_[2])) +
            (4.272*(-3 + complex(0,1)*x_[3] + x_[2]))/(-1 + complex(0,1)*x_[3] + x_[2]) + (x_[1]*(0. + (
                complex(0.06721427746115585,0.0031882562253096425) - complex(1.7472374660574468,-0.008345192421516656)*x_[1] + 0.534*x_[1]**2)/(
                    complex(-0.0014272717944412566,-1.793739297851695e-21) + complex(0.16855723761492827,0.005951877656781462)*x_[1] -
                    complex(3.2849786279185222,-0.015627693218553906)*x_[1]**2 + 1.*x_[1]**3))*(-1 + complex(0,1)*x_[3] + x_[2]))/(complex(0,1)*x_[3] + x_[2]))


class GUTv1_theory(Model):
    def __init__(self, pars):
        super(GUTv1_theory, self).__init__(pars)

    @property
    def val(self):
        """ Implement and return model prediction here using self.x_ as parameter point
        Convention: self.x_ = [mR, r1, Rer2, Imr2]
        """
        return (t11(self.x_), t12(self.x_), t13(self.x_), t22(self.x_), t23(self.x_), t33(self.x_))


@jit
def d11(x_):
    return 0.8251461863096498*(0.8251461863096498*x_[0]) + 0.544910563973106*(0.544910563973106*cmath.sqrt(0.0000742 + x_[0]**2)) - complex(0.14393995936415277,-0.03856859587019336)*(0. - complex(0.14393995936415277,-0.03856859587019336)*cmath.sqrt(0.0025139999999999997 + x_[0]**2))

@jit
def d12(x_):
    return complex(-0.27087999939234897,0.02428963334345038)*(0.8251461863096498*x_[0]) + complex(0.6073208172717925,0.016040403535128984)*(0.544910563973106*cmath.sqrt(0.0000742 + x_[0]**2)) + 0.7462829021249207*(0. - complex(0.14393995936415277,-0.03856859587019336)*cmath.sqrt(0.0025139999999999997 + x_[0]**2))


@jit
def d13(x_):
    return complex(0.49469382668102385,0.021114656131880235)*(0.8251461863096498*x_[0]) - complex(0.5777388515436382,-0.01394371006230815)*(0.544910563973106*cmath.sqrt(0.0000742 + x_[0]**2)) + 0.6487338294761988*(0. - complex(0.14393995936415277,-0.03856859587019336)*cmath.sqrt(0.0025139999999999997 + x_[0]**2))


@jit
def d22(x_):
    return complex(-0.27087999939234897,0.02428963334345038)*( - complex(0.27087999939234897,-0.02428963334345038)*x_[0]) + complex(0.6073208172717925,0.016040403535128984)*(complex(0.6073208172717925,0.016040403535128984)*cmath.sqrt(0.0000742 + x_[0]**2)) + 0.7462829021249207*(0.7462829021249207*cmath.sqrt(0.0025139999999999997 + x_[0]**2))

@jit
def d23(x_):
    return complex(0.49469382668102385,0.021114656131880235)*(- complex(0.27087999939234897,-0.02428963334345038)*x_[0]) - complex(0.5777388515436382,-0.01394371006230815)*(complex(0.6073208172717925,0.016040403535128984)*cmath.sqrt(0.0000742 + x_[0]**2)) + 0.6487338294761988*(0.7462829021249207*cmath.sqrt(0.0025139999999999997 + x_[0]**2))


@jit
def d33(x_):
    return complex(0.49469382668102385,0.021114656131880235)*(complex(0.49469382668102385,0.021114656131880235)*x_[0]) - complex(0.5777388515436382,-0.01394371006230815)*(- complex(0.5777388515436382,-0.01394371006230815)*cmath.sqrt(0.0000742 + x_[0]**2)) + 0.6487338294761988*(0.6487338294761988*cmath.sqrt(0.0025139999999999997 + x_[0]**2))





class GUTv1_data(Model):
    def __init__(self, pars):
        super(GUTv1_data, self).__init__(pars)

    @property
    def val(self):
        """ Implement and return model prediction here using self.x_ as parameter point
        Convention: self.x_ = [m1]
        """
        return (d11(self.x_), d12(self.x_), d13(self.x_), d22(self.x_), d23(self.x_), d33(self.x_))

if __name__ == "__main__":
    m = GUTv1_data( [1] )
    print(m)
    print(m((1)))


    m = GUTv1_theory( [  1, [1,10], [1,3], 2] )
    print(m)
    print(  m( (1,2)) )




