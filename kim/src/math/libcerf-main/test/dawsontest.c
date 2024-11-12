/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File dawsontest.c
 *   Test the complex Dawson function.
 *
 * Copyright:
 *   (C) 2012 Massachusetts Institute of Technology
 *   (C) 2013 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   ../LICENSE
 *
 * Authors:
 *   Steven G. Johnson, Massachusetts Institute of Technology, 2012, core author
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013, package maintainer
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 *
 * Revision history:
 *   ../CHANGELOG
 */

#include "cerf.h"
#include "testtool.h"

int main()
{
    result_t result = {0, 0};
    // dawson(z[i]), evaluated with Maple
    ZTEST(
        result, 1e-13, cdawson(C(2, 1)),
        C(0.1635394094345355614904345232875688576839, -0.1531245755371229803585918112683241066853));
    ZTEST(
        result, 1e-13, cdawson(C(-2, 1)),
        C(-0.1635394094345355614904345232875688576839,
          -0.1531245755371229803585918112683241066853));
    ZTEST(
        result, 1e-13, cdawson(C(2, -1)),
        C(0.1635394094345355614904345232875688576839, 0.1531245755371229803585918112683241066853));
    ZTEST(
        result, 1e-13, cdawson(C(-2, -1)),
        C(-0.1635394094345355614904345232875688576839, 0.1531245755371229803585918112683241066853));
    ZTEST(
        result, 1e-13, cdawson(C(-28, 9)),
        C(-0.01619082256681596362895875232699626384420,
          -0.005210224203359059109181555401330902819419));
    ZTEST(
        result, 1e-13, cdawson(C(33, -21)),
        C(0.01078377080978103125464543240346760257008,
          0.006866888783433775382193630944275682670599));
    ZTEST(
        result, 1e-13, cdawson(C(1e3, 1e3)),
        C(-0.5808616819196736225612296471081337245459, 0.6688593905505562263387760667171706325749));
    ZTEST(result, 1e-13, cdawson(C(-1000, -3001)), C(Inf, -Inf));
    ZTEST(
        result, 1e-13, cdawson(C(1e-8, 5.1e-3)),
        C(0.1000052020902036118082966385855563526705e-7,
          0.005100088434920073153418834680320146441685));
    ZTEST(
        result, 1e-13, cdawson(C(4.95e-3, -4.9e-3)),
        C(0.004950156837581592745389973960217444687524,
          -0.004899838305155226382584756154100963570500));
    ZTEST(
        result, 1e-13, cdawson(C(5.1e-3, 5.1e-3)),
        C(0.005100176864319675957314822982399286703798,
          0.005099823128319785355949825238269336481254));
    ZTEST(
        result, 1e-13, cdawson(C(0.5, 4.9e-3)),
        C(0.4244534840871830045021143490355372016428,
          0.002820278933186814021399602648373095266538));
    ZTEST(
        result, 1e-13, cdawson(C(-0.5e1, 4.9e-4)),
        C(-0.1021340733271046543881236523269967674156,
          -0.00001045696456072005761498961861088944159916));
    ZTEST(
        result, 1e-13, cdawson(C(-0.5e2, -4.9e-5)),
        C(-0.01000200120119206748855061636187197886859,
          0.9805885888237419500266621041508714123763e-8));
    ZTEST(
        result, 1e-13, cdawson(C(0.5e3, 4.9e-6)),
        C(0.001000002000012000023960527532953151819595,
          -0.9800058800588007290937355024646722133204e-11));
    ZTEST(
        result, 1e-13, cdawson(C(0.5, 5.1e-3)),
        C(0.4244549085628511778373438768121222815752,
          0.002935393851311701428647152230552122898291));
    ZTEST(
        result, 1e-13, cdawson(C(-0.5e1, 5.1e-4)),
        C(-0.1021340732357117208743299813648493928105,
          -0.00001088377943049851799938998805451564893540));
    ZTEST(
        result, 1e-13, cdawson(C(-0.5e2, -5.1e-5)),
        C(-0.01000200120119126652710792390331206563616,
          0.1020612612857282306892368985525393707486e-7));
    ZTEST(
        result, 1e-13, cdawson(C(1e-6, 2e-6)),
        C(0.1000000000007333333333344266666666664457e-5,
          0.2000000000001333333333323199999999978819e-5));
    ZTEST(result, 1e-13, cdawson(C(2e-6, 0)), C(0.1999999999994666666666675199999999990248e-5, 0));
    ZTEST(result, 1e-13, cdawson(C(2, 0)), C(0.3013403889237919660346644392864226952119, 0));
    ZTEST(result, 1e-13, cdawson(C(20, 0)), C(0.02503136792640367194699495234782353186858, 0));
    ZTEST(result, 1e-13, cdawson(C(200, 0)), C(0.002500031251171948248596912483183760683918, 0));
    ZTEST(result, 1e-13, cdawson(C(0, 4.9e-3)), C(0, 0.004900078433419939164774792850907128053308));
    ZTEST(
        result, 1e-13, cdawson(C(0, -5.1e-3)), C(0, -0.005100088434920074173454208832365950009419));
    ZTEST(result, 1e-13, cdawson(C(0, 2e-6)), C(0, 0.2000000000005333333333341866666666676419e-5));
    ZTEST(result, 1e-13, cdawson(C(0, -2)), C(0, -48.16001211429122974789822893525016528191));
    ZTEST(result, 1e-13, cdawson(C(0, 20)), C(0, 0.4627407029504443513654142715903005954668e174));
    ZTEST(result, 1e-13, cdawson(C(0, -200)), C(0, -Inf));
    ZTEST(result, 1e-13, cdawson(C(Inf, 0)), C(0, 0));
    ZTEST(result, 1e-13, cdawson(C(-Inf, 0)), C(-0, 0));
    // ZTEST(result, 1e-13, cdawson(C(0, Inf)), C(0, Inf));
    // ZTEST(result, 1e-13, cdawson(C(0, -Inf)), C(0, -Inf));
    ZTEST(result, 1e-13, cdawson(C(Inf, Inf)), C(NaN, NaN));
    ZTEST(result, 1e-13, cdawson(C(Inf, -Inf)), C(NaN, NaN));
    ZTEST(result, 1e-13, cdawson(C(NaN, NaN)), C(NaN, NaN));
    ZTEST(result, 1e-13, cdawson(C(NaN, 0)), C(NaN, 0));
    ZTEST(result, 1e-13, cdawson(C(0, NaN)), C(0, NaN));
    ZTEST(result, 1e-13, cdawson(C(NaN, Inf)), C(NaN, NaN));
    ZTEST(result, 1e-13, cdawson(C(Inf, NaN)), C(NaN, NaN));
    ZTEST(
        result, 1e-13, cdawson(C(39, 6.4e-5)),
        C(0.01282473148489433743567240624939698290584,
          -0.2105957276516618621447832572909153498104e-7));
    ZTEST(
        result, 1e-13, cdawson(C(41, 6.09e-5)),
        C(0.01219875253423634378984109995893708152885,
          -0.1813040560401824664088425926165834355953e-7));
    ZTEST(
        result, 1e-13, cdawson(C(4.9e7, 5e-11)),
        C(0.1020408163265306334945473399689037886997e-7,
          -0.1041232819658476285651490827866174985330e-25));
    ZTEST(
        result, 1e-13, cdawson(C(5.1e7, 4.8e-11)),
        C(0.9803921568627452865036825956835185367356e-8,
          -0.9227220299884665067601095648451913375754e-26));
    ZTEST(
        result, 1e-13, cdawson(C(1e9, 2.4e-12)),
        C(0.5000000000000000002500000000000000003750e-9,
          -0.1200000000000000001800000188712838420241e-29));
    ZTEST(
        result, 1e-13, cdawson(C(1e11, 2.4e-14)),
        C(5.00000000000000000000025000000000000000000003e-12,
          -1.20000000000000000000018000000000000000000004e-36));
    ZTEST(
        result, 1e-13, cdawson(C(1e13, 2.4e-16)),
        C(5.00000000000000000000000002500000000000000000e-14,
          -1.20000000000000000000000001800000000000000000e-42));
    ZTEST(result, 1e-13, cdawson(C(1e300, 2.4e-303)), C(5e-301, 0));

    printf("%i/%i tests failed\n", result.failed, result.total);
    return result.failed;
}
