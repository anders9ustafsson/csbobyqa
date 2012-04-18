using System;
using System.Collections.Generic;

namespace Cureos.Numerics
{
    public static class Bobyqa
    {
        #region FIELDS

        private const double ONEMIN = -1.0;
        private const double ZERO = 0.0;
        private const double TENTH = 0.1;
        private const double HALF = 0.5;
        private const double ONE = 1.0;
        private const double TWO = 2.0;
        private const double TEN = 10.0;

        private static readonly string LF = Environment.NewLine;
        private static readonly string L320 = LF + "Return from BOBYQA because of much cancellation in a denominator.";
        private static readonly string L390 = LF + "Return from BOBYQA because CALFUN has been called MAXFUN times.";
        private static readonly string L400 = LF + "Function number {0:I6}    F ={1,18:E10}" + LF +
                                              "The corresponding X is:{2}";
        private static readonly string L710 = LF + "Least value of F ={0,15:E9}" + LF + "The corresponding X is:{1}";

        #endregion

        #region METHODS

        private static void BOBYQA(Func<int, double[], double> calfun, int n, int npt, double[] x,
                                   double[] xl, double[] xu, double rhobeg, double rhoend, int iprint, int maxfun)
        {
            //     This subroutine seeks the least value of a function of many variables,
            //     by applying a trust region method that forms quadratic models by
            //     interpolation. There is usually some freedom in the interpolation
            //     conditions, which is taken up by minimizing the Frobenius norm of
            //     the change to the second derivative of the model, beginning with the
            //     zero matrix. The values of the variables are constrained by upper and
            //     lower bounds. The arguments of the subroutine are as follows.
            //
            //     N must be set to the number of variables and must be at least two.
            //     NPT is the number of interpolation conditions. Its value must be in
            //       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
            //       recommended.
            //     Initial values of the variables must be set in X(1),X(2),...,X(N). They
            //       will be changed to the values that give the least calculated F.
            //     For I=1,2,...,N, XL[I] and XU[I] must provide the lower and upper
            //       bounds, respectively, on X[I]. The construction of quadratic models
            //       requires XL[I] to be strictly less than XU[I] for each I. Further,
            //       the contribution to a model from changes to the I-th variable is
            //       damaged severely by rounding errors if XU[I]-XL[I] is too small.
            //     RHOBEG and RHOEND must be set to the initial and final values of a trust
            //       region radius, so both must be positive with RHOEND no greater than
            //       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
            //       expected change to a variable, while RHOEND should indicate the
            //       accuracy that is required in the final values of the variables. An
            //       error return occurs if any of the differences XU[I]-XL[I], I=1,...,N,
            //       is less than 2*RHOBEG.
            //     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
            //       amount of printing. Specifically, there is no output if IPRINT=0 and
            //       there is output only at the return if IPRINT=1. Otherwise, each new
            //       value of RHO is printed, with the best vector of variables so far and
            //       the corresponding value of the objective function. Further, each new
            //       value of F with its variables are output if IPRINT=3.
            //     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
            //
            //     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
            //     F to the value of the objective function for the current values of the
            //     variables X(1),X(2),...,X(N), which are generated automatically in a
            //     way that satisfies the bounds given in XL and XU.

            //     Return if the value of NPT is unacceptable.

            var np = n + 1;
            if (npt < n + 2 || npt > ((n + 2) * np) / 2)
            {
                Console.WriteLine(LF + "Return from BOBYQA because NPT is not in the required interval");
                return;
            }

            var ndim = npt + n;

            var sl = new double[1 + n];
            var su = new double[1 + n];

            //     Return if there is insufficient space between the bounds. Modify the
            //     initial X if necessary in order to avoid conflicts between the bounds
            //     and the construction of the first quadratic model. The lower and upper
            //     bounds on moves from the updated X are set now, in the ISL and ISU
            //     partitions of W, in order to provide useful and exact information about
            //     components of X that become within distance RHOBEG from their bounds.

            for (var j = 1; j <= n; ++j)
            {
                double temp = xu[j] - xl[j];
                if (temp < rhobeg + rhobeg)
                {
                    Console.WriteLine(LF + "Return from BOBYQA because one of the differences XU[I]-XL[I] is less than 2*RHOBEG.");
                    return;
                }
                sl[j] = xl[j] - x[j];
                su[j] = xu[j] - x[j];
                if (sl[j] >= -rhobeg)
                {
                    if (sl[j] >= ZERO)
                    {
                        x[j] = xl[j];
                        sl[j] = ZERO;
                        su[j] = temp;
                    }
                    else
                    {
                        x[j] = xl[j] + rhobeg;
                        sl[j] = -rhobeg;
                        su[j] = Math.Max(xu[j] - x[j], rhobeg);
                    }
                }
                else if (su[j] <= rhobeg)
                {
                    if (su[j] <= ZERO)
                    {
                        x[j] = xu[j];
                        sl[j] = -temp;
                        su[j] = ZERO;
                    }
                    else
                    {
                        x[j] = xu[j] - rhobeg;
                        sl[j] = Math.Min(xl[j] - x[j], -rhobeg);
                        su[j] = rhobeg;
                    }
                }
            }

            //     Make the call of BOBYQB.
            BOBYQB(calfun, n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun, ndim, sl, su);
        }

        private static void BOBYQB(Func<int, double[], double> CALFUN, int N, int NPT, double[] X, double[] XL,
                                   double[] XU, double RHOBEG, double RHOEND, int IPRINT, int MAXFUN, int NDIM,
                                   double[] SL, double[] SU)
        {
            //     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN
            //       are identical to the corresponding arguments in SUBROUTINE BOBYQA.
            //     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
            //       All the components of every XOPT are going to satisfy the bounds
            //       SL[I] .LEQ. XOPT[I] .LEQ. SU[I], with appropriate equalities when
            //       XOPT is on a constraint boundary.

            //     Set some constants.

            var NP = N + 1;
            var NPTM = NPT - NP;
            var NH = (N * NP) / 2;

            //     XBASE holds a shift of origin that should reduce the contributions
            //       from rounding errors to values of the model and Lagrange functions.
            //     XPT is a two-dimensional array that holds the coordinates of the
            //       interpolation points relative to XBASE.
            //     FVAL holds the values of F at the interpolation points.
            //     XOPT is set to the displacement from XBASE of the trust region centre.
            //     GOPT holds the gradient of the quadratic model at XBASE+XOPT.
            //     HQ holds the explicit second derivatives of the quadratic model.
            //     PQ contains the parameters of the implicit second derivatives of the
            //       quadratic model.
            //     BMAT holds the last N columns of H.
            //     ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
            //       this factorization being ZMAT times ZMAT^T, which provides both the
            //       correct rank and positive semi-definiteness.
            //     NDIM is the first dimension of BMAT and has the value NPT+N.
            //     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
            //       vector of variables for the next call of CALFUN. XNEW also satisfies
            //       the SL and SU constraints in the way that has just been mentioned.
            //     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
            //       in order to increase the denominator in the updating of UPDATE.
            //     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
            //     VLAG contains the values of the Lagrange functions at a new point X.
            //       They are part of a product that requires VLAG to be of length NDIM.
            //     W is a one-dimensional array that is used for working space. Its length
            //       must be at least 3*NDIM = 3*(NPT+N).

            var XBASE = new double[1 + N];
            var XPT = new double[1 + NPT,1 + N];
            var FVAL = new double[1 + NPT];
            var XOPT = new double[1 + N];
            var GOPT = new double[1 + N];
            var HQ = new double[1 + N * NP / 2];
            var PQ = new double[1 + NPT];
            var BMAT = new double[1 + NDIM,1 + N];
            var ZMAT = new double[1 + NPT,1 + NPT - NP];
            var XNEW = new double[1 + N];
            var XALT = new double[1 + N];
            var D = new double[1 + N];
            var VLAG = new double[1 + NDIM];

            double F, DIFFC, DISTSQ, RATIO;

            //     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
            //     BMAT and ZMAT for the first iteration, with the corresponding values of
            //     of NF and KOPT, which are the number of calls of CALFUN so far and the
            //     index of the interpolation point at the trust region centre. Then the
            //     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
            //     less than NPT. GOPT will be updated if KOPT is different from KBASE.

            var NF = 0;
            var KOPT = 0;
            PRELIM(N, NPT, X, XL, XU, RHOBEG, IPRINT, MAXFUN, XBASE, XPT, FVAL, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL, SU,
                   ref NF, ref KOPT);
            var XOPTSQ = ZERO;
            for (var i = 1; i <= N; ++i)
            {
                XOPT[i] = XPT[KOPT, i];
                XOPTSQ += XOPT[i] * XOPT[i];
            }
            var FSAVE = FVAL[1];
            if (NF < NPT)
            {
                if (IPRINT > 0)
                    Console.WriteLine(L390);
                goto L_720;
            }
            var KBASE = 1;


            //     Complete the settings that are required for the iterative procedure.

            var RHO = RHOBEG;
            var DELTA = RHO;
            var NRESC = NF;
            var NTRITS = 0;
            var DIFFA = ZERO;
            var DIFFB = ZERO;
            var ITEST = 0;
            var NFSAV = NF;

            //     Update GOPT if necessary before the first iteration and after each
            //     call of RESCUE that makes a call of CALFUN.

            L_20:
            if (KOPT != KBASE)
            {
                var ih = 0;
                for (var j = 1; j <= N; ++j)
                {
                    for (var i = 1; i <= j; ++i)
                    {
                        ih = ih + 1;
                        if (i < j) GOPT[j] += HQ[ih] * XOPT[i];
                        GOPT[i] += HQ[ih] * XOPT[j];
                    }
                }
                if (NF > NPT)
                {
                    for (var k = 1; k <= NPT; ++k)
                    {
                        var temp = ZERO;
                        for (var j = 1; j <= N; ++j) temp += XPT[k, j] * XOPT[j];
                        temp *= PQ[k];
                        for (var i = 1; i <= N; ++i) GOPT[i] += temp * XPT[k, i];
                    }
                }
            }

            //     Generate the next point in the trust region that provides a small value
            //     of the quadratic model subject to the constraints on the variables.
            //     The integer NTRITS is set to the number "trust region" iterations that
            //     have occurred since the last "alternative" iteration. If the length
            //     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
            //     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.

            L_60:

            double DSQ, CRVMIN;
            var gnew = new double[1 + N];
            TRSBOX(N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL, SU, DELTA, XNEW, D, gnew, out DSQ, out CRVMIN);
            var DNORM = Math.Min(DELTA, Math.Sqrt(DSQ));
            if (DNORM < HALF * RHO)
            {
                NTRITS = -1;
                DISTSQ = TEN * TEN * RHO * RHO;
                if (NF <= NFSAV + 2) goto L_650;

                //     The following choice between labels 650 and 680 depends on whether or
                //     not our work with the current RHO seems to be complete. Either RHO is
                //     decreased or termination occurs if the errors in the quadratic model at
                //     the last three interpolation points compare favourably with predictions
                //     of likely improvements to the model within distance HALF*RHO of XOPT.
                //
                var errbig = Math.Max(DIFFA, Math.Max(DIFFB, DIFFC));
                var frhosq = 0.125 * RHO * RHO;
                if (CRVMIN > ZERO && errbig > frhosq * CRVMIN) goto L_650;
                var bdtol = errbig / RHO;
                for (var j = 1; j <= N; ++j)
                {
                    var bdtest = bdtol;
                    if (XNEW[j] == SL[j]) bdtest = gnew[j];
                    if (XNEW[j] == SU[j]) bdtest = -gnew[j];
                    if (bdtest < bdtol)
                    {
                        var CURV = HQ[(j + j * j) / 2];
                        for (var k = 1; k <= NPT; ++k)
                        {
                            CURV = CURV + PQ[k] * XPT[k, j] * XPT[k, j];
                        }
                        bdtest = bdtest + HALF * CURV * RHO;
                        if (bdtest < bdtol) goto L_650;
                    }
                }
                goto L_680;
            }
            ++NTRITS;

            //     Severe cancellation is likely to occur if XOPT is too far from XBASE.
            //     If the following test holds, then XBASE is shifted so that XOPT becomes
            //     zero. The appropriate changes are made to BMAT and to the second
            //     derivatives of the current model, beginning with the changes to BMAT
            //     that do not depend on ZMAT. VLAG is used temporarily for working space.

            L_90:

            if (DSQ <= 1.0E-3 * XOPTSQ)
            {
                var nWork = new double[1 + N];
                var nptWork = new double[1 + NPT];

                var fracsq = 0.25 * XOPTSQ;
                var sumpq = ZERO;
                for (var k = 1; k <= NPT; ++k)
                {
                    sumpq += PQ[k];
                    var sum = -HALF * XOPTSQ;
                    for (var i = 1; i <= N; ++i) sum += XPT[k, i] * XOPT[i];
                    nptWork[k] = sum;
                    var temp = fracsq - HALF * sum;
                    for (var i = 1; i <= N; ++i)
                    {
                        nWork[i] = BMAT[k, i];
                        VLAG[i] = sum * XPT[k, i] + temp * XOPT[i];
                        var ip = NPT + i;
                        for (var j = 1; j <= i; ++j) BMAT[ip, j] += nWork[i] * VLAG[j] + VLAG[i] * nWork[j];
                    }
                }

                //     Then the revisions of BMAT that depend on ZMAT are calculated.

                for (var jj = 1; jj <= NPTM; ++jj)
                {
                    var sumz = ZERO;
                    var sumw = ZERO;
                    for (var k = 1; k <= NPT; ++k)
                    {
                        sumz += ZMAT[k, jj];
                        VLAG[k] = nptWork[k] * ZMAT[k, jj];
                        sumw += VLAG[k];
                    }
                    for (var j = 1; j <= N; ++j)
                    {
                        var sum = (fracsq * sumz - HALF * sumw) * XOPT[j];
                        for (var k = 1; k <= NPT; ++k)
                        {
                            sum += VLAG[k] * XPT[k, j];
                        }
                        nWork[j] = sum;
                        for (var k = 1; k <= NPT; ++k) BMAT[k, j] += sum * ZMAT[k, jj];
                    }
                    for (var i = 1; i <= N; ++i)
                    {
                        var ip = i + NPT;
                        var temp = nWork[i];
                        for (var j = 1; j <= i; ++j) BMAT[ip, j] += temp * nWork[j];
                    }
                }

                //     The following instructions complete the shift, including the changes
                //     to the second derivative parameters of the quadratic model.

                var ih = 0;
                for (var j = 1; j <= N; ++j)
                {
                    nWork[j] = -HALF * sumpq * XOPT[j];
                    for (var k = 1; k <= NPT; ++k)
                    {
                        nWork[j] += PQ[k] * XPT[k, j];
                        XPT[k, j] -= XOPT[j];
                    }
                    for (var i = 1; i <= j; ++i)
                    {
                        ++ih;
                        HQ[ih] += nWork[i] * XOPT[j] + XOPT[i] * nWork[j];
                        BMAT[NPT + i, j] = BMAT[NPT + j, i];
                    }
                }
                for (var i = 1; i <= N; ++i)
                {
                    XBASE[i] += XOPT[i];
                    XNEW[i] -= XOPT[i];
                    SL[i] -= XOPT[i];
                    SU[i] -= XOPT[i];
                    XOPT[i] = ZERO;
                }
                XOPTSQ = ZERO;
            }

            if (NTRITS == 0) goto L_210;
            goto L_230;

            //     XBASE is also moved to XOPT by a call of RESCUE. This calculation is
            //     more expensive than the previous shift, because new matrices BMAT and
            //     ZMAT are generated from scratch, which may include the replacement of
            //     interpolation points whose positions seem to be causing near linear
            //     dependence in the interpolation conditions. Therefore RESCUE is called
            //     only if rounding errors have reduced by at least a factor of two the
            //     denominator of the formula for updating the H matrix. It provides a
            //     useful safeguard, but is not invoked in most applications of BOBYQA.

            L_190:

            NFSAV = NF;
            KBASE = KOPT;
            RESCUE(N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL,
                   XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL, SU, ref NF, DELTA, ref KOPT,
                   VLAG);

            //     XOPT is updated now in case the branch below to label 720 is taken.
            //     Any updating of GOPT occurs after the branch below to label 20, which
            //     leads to a trust region iteration as does the branch to label 60.

            XOPTSQ = ZERO;
            if (KOPT != KBASE)
            {
                for (var I = 1; I <= N; ++I)
                {
                    XOPT[I] = XPT[KOPT, I];
                    XOPTSQ = XOPTSQ + XOPT[I] * XOPT[I];
                }
            }
            if (NF < 0)
            {
                NF = MAXFUN;
                if (IPRINT > 0) Console.WriteLine(L390);
                goto L_720;
            }
            NRESC = NF;
            if (NFSAV < NF)
            {
                NFSAV = NF;
                goto L_20;
            }
            if (NTRITS > 0) goto L_60;

            //     Pick two alternative vectors of variables, relative to XBASE, that
            //     are suitable as new positions of the KNEW-th interpolation point.
            //     Firstly, XNEW is set to the point on a line through XOPT and another
            //     interpolation point that minimizes the predicted value of the next
            //     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
            //     and SU bounds. Secondly, XALT is set to the best feasible point on
            //     a constrained version of the Cauchy step of the KNEW-th Lagrange
            //     function, the corresponding value of the square of this function
            //     being returned in CAUCHY. The choice between these alternatives is
            //     going to be made when the denominator is calculated.

            L_210:

            int KNEW;
            double ADELT, ALPHA, CAUCHY;

            ALTMOV(N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL, SU, KOPT,
                   KNEW, ADELT, XNEW, XALT, out ALPHA, out CAUCHY);
            for (var i = 1; i <= N; ++i) D[i] = XNEW[i] - XOPT[i];

            //     Calculate VLAG and BETA for the current choice of D. The scalar
            //     product of D with XPT(K,.) is going to be held in W(NPT+K) for
            //     use when VQUAD is calculated.

            L_230:

            var W = new double[1 + 2 * NPT];

            for (var k = 1; k <= NPT; ++k)
            {
                var suma = ZERO;
                var sumb = ZERO;
                var sum = ZERO;
                for (var j = 1; j <= N; ++j)
                {
                    suma += XPT[k, j] * D[j];
                    sumb += XPT[k, j] * XOPT[j];
                    sum += BMAT[k, j] * D[j];
                }
                W[k] = suma * (HALF * suma + sumb);
                VLAG[k] = sum;
                W[NPT + k] = suma;
            }
            var BETA = ZERO;
            for (var jj = 1; jj >= NPTM; ++jj)
            {
                var sum = ZERO;
                for (var k = 1; k <= NPT; ++k) sum += ZMAT[k, jj] * W[k];
                BETA = BETA - sum * sum;
                for (var k = 1; k <= NPT; ++k) VLAG[k] = VLAG[k] + sum * ZMAT[k, jj];
            }
            DSQ = ZERO;
            var bsum = ZERO;
            var dx = ZERO;
            for (var j = 1; j <= N; ++j)
            {
                DSQ = DSQ + D[j] * D[j];
                var sum = ZERO;
                for (var k = 1; k <= NPT; ++k) sum += W[k] * BMAT[k, j];
                bsum += sum * D[j];
                var jp = NPT + j;
                for (var i = 1; i <= N; ++i) sum += BMAT[jp, i] * D[i];
                VLAG[jp] = sum;
                bsum += sum * D[j];
                dx += D[j] * XOPT[j];
            }
            BETA = dx * dx + DSQ * (XOPTSQ + dx + dx + HALF * DSQ) + BETA - bsum;
            VLAG[KOPT] += ONE;

            //     If NTRITS is zero, the denominator may be increased by replacing
            //     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
            //     rounding errors have damaged the chosen denominator.

            double DENOM;
            if (NTRITS == 0)
            {
                DENOM = VLAG[KNEW] * VLAG[KNEW] + ALPHA * BETA;
                if (DENOM < CAUCHY && CAUCHY > ZERO)
                {
                    for (var i = 1; i <= N; ++i)
                    {
                        XNEW[i] = XALT[i];
                        D[i] = XNEW[i] - XOPT[i];
                    }
                    CAUCHY = ZERO;
                    goto L_230;
                }
                if (DENOM <= HALF * VLAG[KNEW] * VLAG[KNEW])
                {
                    if (NF > NRESC) goto L_190;
                    if (IPRINT > 0)
                        Console.WriteLine(L320);
                    goto L_720;
                }
            }

                //     Alternatively, if NTRITS is positive, then set KNEW to the index of
                //     the next interpolation point to be deleted to make room for a trust
                //     region step. Again RESCUE may be called if rounding errors have damaged
                //     the chosen denominator, which is the reason for attempting to select
                //     KNEW before calculating the next value of the objective function.

            else
            {
                var DELSQ = DELTA * DELTA;
                var SCADEN = ZERO;
                var BIGLSQ = ZERO;
                KNEW = 0;
                for (var k = 1; k <= NPT; ++k)
                {
                    if (k == KOPT) continue;
                    var HDIAG = ZERO;
                    for (var jj = 1; jj <= NPTM; ++jj) HDIAG = HDIAG + ZMAT[k, jj] * ZMAT[k, jj];
                    var DEN = BETA * HDIAG + VLAG[k] * VLAG[k];
                    DISTSQ = ZERO;
                    for (var j = 1; j <= N; ++j) DISTSQ += Math.Pow(XPT[k, j] - XOPT[j], 2.0);
                    var temp = Math.Max(ONE, Math.Pow(DISTSQ / DELSQ, 2.0));
                    if (temp * DEN > SCADEN)
                    {
                        SCADEN = temp * DEN;
                        KNEW = k;
                        DENOM = DEN;
                    }
                    BIGLSQ = Math.Max(BIGLSQ, temp * VLAG[k] * VLAG[k]);
                }
                if (SCADEN <= HALF * BIGLSQ)
                {
                    if (NF > NRESC) goto L_190;
                    if (IPRINT > 0) Console.WriteLine(L320);
                    goto L_720;
                }
            }

            //     Put the variables for the next calculation of the objective function
            //       in XNEW, with any adjustments for the bounds.

            //     Calculate the value of the objective function at XBASE+XNEW, unless
            //       the limit on the number of calculations of F has been reached.

            L_360:

            for (var i = 1; i <= N; ++i)
            {
                X[i] = Math.Min(Math.Max(XL[i], XBASE[i] + XNEW[i]), XU[i]);
                if (XNEW[i] == SL[i]) X[i] = XL[i];
                if (XNEW[i] == SU[i]) X[i] = XU[i];
            }
            if (NF >= MAXFUN)
            {
                if (IPRINT > 0) Console.WriteLine(L390);
                goto L_720;
            }

            ++NF;
            F = CALFUN(N, X);

            if (IPRINT == 3) Console.WriteLine(L400, NF, F, X.PART(1, N).FORMAT());
            if (NTRITS == -1)
            {
                FSAVE = F;
                goto L_720;
            }

            //     Use the quadratic model to predict the change in F due to the step D,
            //       and set DIFF to the error of this prediction.

            var FOPT = FVAL[KOPT];
            var VQUAD = ZERO;
            {
                var ih = 0;
                for (var j = 1; j <= N; ++j)
                {
                    VQUAD += D[j] * GOPT[j];
                    for (var i = 1; i <= j; ++i)
                        VQUAD += HQ[++ih] * (i == j ? HALF : ONE) * D[i] * D[j];
                }
            }
            for (var k = 1; k <= NPT; ++k) VQUAD += HALF * PQ[k] * W[NPT + k] * W[NPT + k];
            var DIFF = F - FOPT - VQUAD;
            DIFFC = DIFFB;
            DIFFB = DIFFA;
            DIFFA = Math.Abs(DIFF);
            if (DNORM > RHO) NFSAV = NF;

            //     Pick the next value of DELTA after a trust region step.

            if (NTRITS > 0)
            {
                if (VQUAD >= ZERO)
                {
                    if (IPRINT > 0)
                        Console.WriteLine(LF + "Return from BOBYQA because a trust region step has failed to reduce Q.");
                    goto L_720;
                }
                RATIO = (F - FOPT) / VQUAD;
                if (RATIO <= TENTH)
                {
                    DELTA = Math.Min(HALF * DELTA, DNORM);
                }
                else if (RATIO <= 0.7)
                {
                    DELTA = Math.Max(HALF * DELTA, DNORM);
                }
                else
                {
                    DELTA = Math.Max(HALF * DELTA, DNORM + DNORM)
                }
                if (DELTA <= 1.5 * RHO) DELTA = RHO;

                //     Recalculate KNEW and DENOM if the new F is less than FOPT.

                if (F < FOPT)
                {
                    var KSAV = KNEW;
                    var DENSAV = DENOM;
                    var DELSQ = DELTA * DELTA;
                    var SCADEN = ZERO;
                    var BIGLSQ = ZERO;
                    KNEW = 0;
                    for (var K = 1; K <= NPT; ++K)
                    {
                        var HDIAG = ZERO;
                        for (var JJ = 1; JJ <= NPTM; ++JJ) HDIAG += ZMAT[K, JJ] * ZMAT[K, JJ];
                        var DEN = BETA * HDIAG + VLAG[K] * VLAG[K];
                        DISTSQ = ZERO;
                        for (var J = 1; J <= N; ++J) DISTSQ += Math.Pow(XPT[K, J] - XNEW[J], 2.0);
                        var TEMP = Math.Max(ONE, Math.Pow(DISTSQ / DELSQ, 2.0));
                        if (TEMP * DEN > SCADEN)
                        {
                            SCADEN = TEMP * DEN;
                            KNEW = K;
                            DENOM = DEN;
                        }
                        BIGLSQ = Math.Max(BIGLSQ, TEMP * VLAG[K] * VLAG[K]);
                    }
                    if (SCADEN <= HALF * BIGLSQ)
                    {
                        KNEW = KSAV;
                        DENOM = DENSAV;
                    }
                }
            }

            //     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
            //     moved. Also update the second derivative terms of the model.

            UPDATE(N, NPT, BMAT, ZMAT, NDIM, VLAG, BETA, DENOM, KNEW);
            var pqold = PQ[KNEW];
            PQ[KNEW] = ZERO;
            {
                var ih = 0;
                for (var i = 1; i <= N; ++i)
                {
                    var temp = pqold * XPT[KNEW, i];
                    for (var j = 1; j <= i; ++j) HQ[++ih] += temp * XPT[KNEW, j];
                }
            }
            for (var jj = 1; jj <= NPTM; ++jj)
            {
                var temp = DIFF * ZMAT[KNEW, jj];
                for (var k = 1; k <= NPT; ++k) PQ[k] += temp * ZMAT[k, jj];
            }

            //     Include the new interpolation point, and make the changes to GOPT at
            //     the old XOPT that are caused by the updating of the quadratic model.

            FVAL[KNEW] = F;
            for (var i = 1; i <= N; ++i)
            {
                XPT[KNEW, i] = XNEW[i];
                W[i] = BMAT[KNEW, i];
            }
            for (var k = 1; k <= NPT; ++k)
            {
                var suma = ZERO;
                for (var jj = 1; jj <= NPTM; ++jj) suma += ZMAT[KNEW, jj] * ZMAT[k, jj];
                var sumb = ZERO;
                for (var j = 1; j <= N; ++j) sumb += XPT[k, j] * XOPT[j];
                var temp = suma * sumb;
                for (var i = 1; i <= N; ++i) W[i] += temp * XPT[k, i];
            }
            for (var i = 1; i <= N; ++i) GOPT[i] += DIFF * W[i];

            //     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.

            if (F < FOPT)
            {
                KOPT = KNEW;
                XOPTSQ = ZERO;
                var ih = 0;
                for (var j = 1; j <= N; ++j)
                {
                    XOPT[j] = XNEW[j];
                    XOPTSQ += XOPT[j] * XOPT[j];
                    for (var i = 1; i <= j; ++i)
                    {
                        ++ih;
                        if (i < j) GOPT[j] += +HQ[ih] * D[i];
                        GOPT[i] += HQ[ih] * D[j];
                    }
                }
                for (var k = 1; k <= NPT; ++k)
                {
                    var temp = ZERO;
                    for (var j = 1; j <= N; ++j) temp += XPT[k, j] * D[j];
                    temp *= PQ[k];
                    for (var i = 1; i <= N; ++i) GOPT[i] += temp * XPT[k, i];
                }
            }

            //     Calculate the parameters of the least Frobenius norm interpolant to
            //     the current data, the gradient of this interpolant at XOPT being put
            //     into VLAG(NPT+I), I=1,2,...,N.

            if (NTRITS > 0)
            {
                for (var k = 1; k <= NPT; ++k)
                {
                    VLAG[k] = FVAL[k] - FVAL[KOPT];
                    W[k] = ZERO;
                }
                for (var J = 1; J <= NPTM; ++J)
                {
                    var sum = ZERO;
                    for (var k = 1; k <= NPT; ++k) sum += ZMAT[k, J] * VLAG[k];
                    for (var k = 1; k <= NPT; ++k) W[k] = W[k] + sum * ZMAT[k, J];
                }
                for (var k = 1; k <= NPT; ++k)
                {
                    var sum = ZERO;
                    for (var j = 1; j <= N; ++j) sum += XPT[k, j] * XOPT[j];
                    W[k + NPT] = W[k];
                    W[k] *= sum;
                }
                var GQSQ = ZERO;
                var GISQ = ZERO;
                for (var i = 1; i <= N; ++i)
                {
                    var SUM = ZERO;
                    for (var k = 1; k <= NPT; ++k) SUM += BMAT[k, i] * VLAG[k] + XPT[k, i] * W[k];
                    if (XOPT[i] == SL[i])
                    {
                        GQSQ += Math.Pow(Math.Min(ZERO, GOPT[i]), 2.0);
                        GISQ += Math.Pow(Math.Min(ZERO, SUM), 2.0);
                    }
                    else if (XOPT[i] == SU[i])
                    {
                        GQSQ += Math.Pow(Math.Max(ZERO, GOPT[i]), 2.0);
                        GISQ += Math.Pow(Math.Max(ZERO, SUM), 2.0);
                    }
                    else
                    {
                        GQSQ += GOPT[i] * GOPT[i];
                        GISQ += SUM * SUM;
                    }
                    VLAG[NPT + i] = SUM;
                }

                //     Test whether to replace the new quadratic model by the least Frobenius
                //     norm interpolant, making the replacement if the test is satisfied.

                ++ITEST;
                if (GQSQ < TEN * GISQ) ITEST = 0;
                if (ITEST >= 3)
                {
                    for (var i = 1; i <= Math.Max(NPT, NH); ++i)
                    {
                        if (i <= N) GOPT[i] = VLAG[NPT + i];
                        if (i <= NPT) PQ[i] = W[NPT + i];
                        if (i <= NH) HQ[i] = ZERO;
                        ITEST = 0;
                    }
                }
            }

//     If a trust region step has provided a sufficient decrease in F, then
//     branch for another trust region calculation. The case NTRITS=0 occurs
//     when the new interpolation point was reached by an alternative step.

            if (NTRITS == 0 || F <= FOPT + TENTH * VQUAD) goto L_60;

//     Alternatively, find out if the interpolation points are close enough
//       to the best point so far.

            DISTSQ = Math.Max(TWO * TWO * DELTA * DELTA, TEN * TEN * RHO * RHO);

            L_650:

            KNEW = 0;
            for (var k = 1; k <= NPT; ++k)
            {
                var SUM = ZERO;
                for (var j = 1; j <= N; ++j) SUM += Math.Pow(XPT[k, j] - XOPT[j], 2.0);
                if (SUM > DISTSQ)
                {
                    KNEW = k;
                    DISTSQ = SUM;
                }
            }

            //     If KNEW is positive, then ALTMOV finds alternative new positions for
            //     the KNEW-th interpolation point within distance ADELT of XOPT. It is
            //     reached via label 90. Otherwise, there is a branch to label 60 for
            //     another trust region iteration, unless the calculations with the
            //     current RHO are complete.

            if (KNEW > 0)
            {
                var dist = Math.Sqrt(DISTSQ);
                if (NTRITS == -1)
                {
                    DELTA = Math.Min(TENTH * DELTA, HALF * dist);
                    if (DELTA <= 1.5 * RHO) DELTA = RHO;
                }
                NTRITS = 0;
                ADELT = Math.Max(Math.Min(TENTH * dist, DELTA), RHO);
                DSQ = ADELT * ADELT;
                goto L_90;
            }
            if (NTRITS == -1) goto L_680;
            if (RATIO > ZERO || Math.Max(DELTA, DNORM) > RHO) goto L_60;
//
//     The calculations with the current value of RHO are complete. Pick the
//       next values of RHO and DELTA.

            L_680:

            if (RHO > RHOEND)
            {
                DELTA = HALF * RHO;
                RATIO = RHO / RHOEND;

                if (RATIO <= 16.0)
                    RHO = RHOEND;
                else if (RATIO <= 250.0)
                    RHO = Math.Sqrt(RATIO) * RHOEND;
                else
                    RHO = TENTH * RHO;

                DELTA = Math.Max(DELTA, RHO);
                if (IPRINT >= 2)
                {
                    var bestX = new double[N];
                    for (var i = 1; i <= N; ++i) bestX[i] = XBASE[i] + XOPT[i];

                    if (IPRINT >= 3) Console.WriteLine();
                    Console.WriteLine("New RHO ={0,11:E4}" + LF + "Number of function values ={1:I6}", RHO, NF);
                    Console.WriteLine(LF + "Least value of F ={0,23:F15}" + LF + "The corresponding X is: {1}",
                                      FVAL[KOPT], bestX.PART(1, N).FORMAT());
                }
                NTRITS = 0;
                NFSAV = NF;
                goto L_60;
            }

            //     Return from the calculation, after another Newton-Raphson step, if
            //       it is too short to have been tried before.

            if (NTRITS == -1) goto L_360;

            L_720:

            if (FVAL[KOPT] <= FSAVE)
            {
                for (var i = 1; i <= N; ++i)
                {
                    X[i] = Math.Min(Math.Max(XL[i], XBASE[i] + XOPT[i]), XU[i]);
                    if (XOPT[i] == SL[i]) X[i] = XL[i];
                    if (XOPT[i] == SU[i]) X[i] = XU[i];
                }
                F = FVAL[KOPT];
            }
            if (IPRINT >= 1)
            {
                Console.WriteLine(LF + "At the return from BOBYQA Number of function values = {0}", NF);
                Console.WriteLine(L710, F, X.PART(1, N).FORMAT());
            }
        }

        private static void ALTMOV(int N, int NPT, double[,] XPT, double[] XOPT, double[,] BMAT, 
            double[,] ZMAT, int NDIM, double[] SL, double[] SU, int KOPT, int KNEW, double ADELT, 
            double[] XNEW, double[] XALT, out double ALPHA, out double CAUCHY)
        {
            //     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
            //       the same meanings as the corresponding arguments of BOBYQB.
            //     KOPT is the index of the optimal interpolation point.
            //     KNEW is the index of the interpolation point that is going to be moved.
            //     ADELT is the current trust region bound.
            //     XNEW will be set to a suitable new position for the interpolation point
            //       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
            //       bounds and it should provide a large denominator in the next call of
            //       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
            //       straight lines through XOPT and another interpolation point.
            //     XALT also provides a large value of the modulus of the KNEW-th Lagrange
            //       function subject to the constraints that have been mentioned, its main
            //       difference from XNEW being that XALT-XOPT is a constrained version of
            //       the Cauchy step within the trust region. An exception is that XALT is
            //       not calculated if all components of GLAG (see below) are zero.
            //     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
            //     CAUCHY will be set to the square of the KNEW-th Lagrange function at
            //       the step XALT-XOPT from XOPT for the vector XALT that is returned,
            //       except that CAUCHY is set to zero if XALT is not calculated.

            //     GLAG is a working space vector of length N for the gradient of the
            //       KNEW-th Lagrange function at XOPT.
            //     HCOL is a working space vector of length NPT for the second derivative
            //       coefficients of the KNEW-th Lagrange function.
            //     W is a working space vector of length 2N that is going to hold the
            //       constrained Cauchy step from XOPT of the Lagrange function, followed
            //       by the downhill version of XALT when the uphill step is calculated.

            var GLAG = new double[1 + N];
            var HCOL = new double[1 + NPT];
            var work3 = new double[1 + 2 * N];

            //     Set the first NPT components of W to the leading elements of the
            //     KNEW-th column of the H matrix.

            var CONST = ONE + Math.Sqrt(2.0);

            for (var K = 1; K <= NPT; ++K) HCOL[K] = ZERO;
            for (var J = 1; J >= NPT - N - 1; ++J)
            {
                var TEMP = ZMAT[KNEW, J];
                for (var K = 1; K <= NPT; ++K) HCOL[K] += TEMP * ZMAT[K, J];
            }
            ALPHA = HCOL[KNEW];
            var HA = HALF * ALPHA;

            // TODO Continue implementation!!!
        }

        private static void PRELIM(int N, int NPT, double[] X, double[] XL, double[] XU,
                                   double RHOBEG, int IPRINT, int MAXFUN, double[] XBASE, double[,] XPT, double[] FVAL,
                                   double[] GOPT, double[] HQ, double[] PQ, double[,] BMAT, double[,] ZMAT,
                                   int NDIM, double[] SL, double[] SU, ref int NF, ref int KOPT)
        {
            //     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
            //       same as the corresponding arguments in SUBROUTINE BOBYQA.
            //     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
            //       are the same as the corresponding arguments in BOBYQB, the elements
            //       of SL and SU being set in BOBYQA.
            //     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
            //       it is set by PRELIM to the gradient of the quadratic model at XBASE.
            //       If XOPT is nonzero, BOBYQB will change it to its usual value later.
            //     NF is maintaned as the number of calls of CALFUN so far.
            //     KOPT will be such that the least calculated value of F so far is at
            //       the point XPT(KOPT,.)+XBASE in the space of the variables.
            //
            //     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
            //     BMAT and ZMAT for the first iteration, and it maintains the values of
            //     NF and KOPT. The vector X is also changed by PRELIM.

            //     Set some constants.

            var RHOSQ = RHOBEG * RHOBEG;
            var RECIP = ONE / RHOSQ;
            var NP = N + 1;

            // TODO Continue implementation!!!
        }

        public static void RESCUE(int N, int NPT, double[] XL, double[] XU, int IPRINT, 
            int MAXFUN, double[] XBASE, double[,] XPT, double[] FVAL, double[] XOPT, double[] GOPT, 
            double[] HQ, double[] PQ, double[,] BMAT, double[,] ZMAT, int NDIM, double[] SL, double[] SU, 
            ref int NF, double DELTA, ref int KOPT, double[] VLAG)
        {

            //     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
            //       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
            //       the corresponding arguments of BOBYQB on the entry to RESCUE.
            //     NF is maintained as the number of calls of CALFUN so far, except that
            //       NF is set to -1 if the value of MAXFUN prevents further progress.
            //     KOPT is maintained so that FVAL(KOPT) is the least calculated function
            //       value. Its correct value must be given on entry. It is updated if a
            //       new least function value is found, but the corresponding changes to
            //       XOPT and GOPT have to be made later by the calling program.
            //     DELTA is the current trust region radius.
            //     VLAG is a working space vector that will be used for the values of the
            //       provisional Lagrange functions at each of the interpolation points.
            //       They are part of a product that requires VLAG to be of length NDIM.
            //     The final elements of BMAT and ZMAT are set in a well-conditioned way
            //       to the values that are appropriate for the new interpolation points.
            //     The elements of GOPT, HQ and PQ are also revised to the values that are
            //       appropriate to the final quadratic model.

            //     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and
            //       PTSAUX(2,J) specify the two positions of provisional interpolation
            //       points when a nonzero step is taken along e_J (the J-th coordinate
            //       direction) through XBASE+XOPT, as specified below. Usually these
            //       steps have length DELTA, but other lengths are chosen if necessary
            //       in order to satisfy the given bounds on the variables.
            //     PTSID is also a working space array. It has NPT components that denote
            //       provisional new positions of the original interpolation points, in
            //       case changes are needed to restore the linear independence of the
            //       interpolation conditions. The K-th point is a candidate for change
            //       if and only if PTSID[K] is nonzero. In this case let p and q be the
            //       integer parts of PTSID[K] and (PTSID[K]-p) multiplied by N+1. If p
            //       and q are both positive, the step from XBASE+XOPT to the new K-th
            //       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
            //       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
            //       p=0, respectively.
            //     The first NDIM+NPT elements of the array W are used for working space. 

            var PTSAUX = new double[1 + 2, 1 + N];
            var PTSID = new double[1 + NPT];
            var W = new double[1 + NDIM + NPT];

            //     Set some constants.
            //
            var NP = N + 1;
            var SFRAC = HALF / NP;
            var NPTM = NPT - NP;
        }

        private static void TRSBOX(int N, int NPT, double[,] XPT, double[] XOPT, double[] GOPT, 
            double[] HQ, double[] PQ, double[] SL, double[] SU, double DELTA,
            double[] XNEW, double[] D, double[] GNEW, out double DSQ, out double CRVMIN)
        {
            //     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
            //       meanings as the corresponding arguments of BOBYQB.
            //     DELTA is the trust region radius for the present calculation, which
            //       seeks a small value of the quadratic model within distance DELTA of
            //       XOPT subject to the bounds on the variables.
            //     XNEW will be set to a new vector of variables that is approximately
            //       the one that minimizes the quadratic model within the trust region
            //       subject to the SL and SU constraints on the variables. It satisfies
            //       as equations the bounds that become active during the calculation.
            //     D is the calculated trial step from XOPT, generated iteratively from an
            //       initial value of zero. Thus XNEW is XOPT+D after the final iteration.
            //     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
            //       when D is updated.
            //     DSQ will be set to the square of the length of XNEW-XOPT.
            //     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise
            //       it is set to the least curvature of H that occurs in the conjugate
            //       gradient searches that are not restricted by any constraints. The
            //       value CRVMIN=-1.0D0 is set, however, if all of these searches are
            //       constrained.
            //
            //     A version of the truncated conjugate gradient is applied. If a line
            //     search is restricted by a constraint, then the procedure is restarted,
            //     the values of the variables that are at their bounds being fixed. If
            //     the trust region boundary is reached, then further changes may be made
            //     to D, each one being in the two dimensional space that is spanned
            //     by the current D and the gradient of Q at XOPT+D, staying on the trust
            //     region boundary. Termination occurs when the reduction in Q seems to
            //     be close to the greatest reduction that can be achieved.

            //     The sign of GOPT[I] gives the sign of the change to the I-th variable
            //     that will reduce Q from its value at XOPT. Thus XBDI[I] shows whether
            //     or not to fix the I-th variable at one of its bounds initially, with
            //     NACT being set to the number of fixed variables. D and GNEW are also
            //     set for the first iteration. DELSQ is the upper bound on the sum of
            //     squares of the free variables. QRED is the reduction in Q so far.

            //     XBDI is a working space vector. For I=1,2,...,N, the element XBDI[I] is
            //       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
            //       I-th variable has become fixed at a bound, the bound being SL[I] or
            //       SU[I] in the case XBDI[I]=-1.0 or XBDI[I]=1.0, respectively. This
            //       information is accumulated during the construction of XNEW.
            //     The arrays S, HS and HRED are also used for working space. They hold the
            //       current search direction, and the changes in the gradient of Q along S
            //       and the reduced D, respectively, where the reduced D is the same as D,
            //       except that the components of the fixed variables are zero.

            var XBDI = new double[1 + N];
            var S = new double[1 + N];
            var HS = new double[1 + N];
            var HRED = new double[1 + N];

            var ITERC = 0;
            var NACT = 0;
            var SQSTP = ZERO;
        }

        private static void UPDATE(int N, int NPT, double[,] BMAT, double[,] ZMAT, int NDIM, double[] VLAG, double BETA, double DENOM, double KNEW)
        {
            //
            //     The arrays BMAT and ZMAT are updated, as required by the new position
            //     of the interpolation point that has the index KNEW. The vector VLAG has
            //     N+NPT components, set on entry to the first NPT and last N components
            //     of the product Hw in equation (4.11) of the Powell (2006) paper on
            //     NEWUOA. Further, BETA is set on entry to the value of the parameter
            //     with that name, and DENOM is set to the denominator of the updating
            //     formula. Elements of ZMAT may be treated as zero if their moduli are
            //     at most ZTEST.

            //     The first NDIM elements of W are used for working space.

            var W = new double[1 + NDIM];

            //     Set some constants.

            var NPTM = NPT - N - 1;
            var ZTEST = ZERO;
            for (var K = 1; K <= NPT; ++K)
                for (var J = 1; J <= NPTM; ++J)
                    ZTEST = Math.Max(ZTEST, Math.Abs(ZMAT[K, J]));
            ZTEST = 1.0E-20 * ZTEST;

            // TODO Continue implementation!!!
        }

        #endregion

        #region PRIVATE SUPPORT METHODS

        private static T[] ROW<T>(this T[,] src, int rowidx)
        {
            var cols = src.GetLength(1);
            var dest = new T[cols];
            for (var col = 0; col < cols; ++col) dest[col] = src[rowidx, col];
            return dest;
        }

        private static T[] COL<T>(this T[,] src, int colidx)
        {
            var rows = src.GetLength(0);
            var dest = new T[rows];
            for (var row = 0; row < rows; ++row) dest[row] = src[row, colidx];
            return dest;
        }

        private static T[] PART<T>(this IList<T> src, int from, int to)
        {
            var dest = new T[to - from + 1];
            var destidx = 0;
            for (var srcidx = from; srcidx <= to; ++srcidx, ++destidx) dest[destidx] = src[srcidx];
            return dest;
        }

        private static string FORMAT(this double[] x)
        {
            return String.Concat(Array.ConvertAll(x, val => String.Format("{0,13:F6}", val)));
        }

        #endregion
    }
}
/*
      SUBROUTINE BOBYQA (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,
     1  MAXFUN,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),XL(*),XU(*),W(*)
//
//     This subroutine seeks the least value of a function of many variables,
//     by applying a trust region method that forms quadratic models by
//     interpolation. There is usually some freedom in the interpolation
//     conditions, which is taken up by minimizing the Frobenius norm of
//     the change to the second derivative of the model, beginning with the
//     zero matrix. The values of the variables are constrained by upper and
//     lower bounds. The arguments of the subroutine are as follows.
//
//     N must be set to the number of variables and must be at least two.
//     NPT is the number of interpolation conditions. Its value must be in
//       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
//       recommended.
//     Initial values of the variables must be set in X(1),X(2),...,X(N). They
//       will be changed to the values that give the least calculated F.
//     For I=1,2,...,N, XL[I] and XU[I] must provide the lower and upper
//       bounds, respectively, on X[I]. The construction of quadratic models
//       requires XL[I] to be strictly less than XU[I] for each I. Further,
//       the contribution to a model from changes to the I-th variable is
//       damaged severely by rounding errors if XU[I]-XL[I] is too small.
//     RHOBEG and RHOEND must be set to the initial and final values of a trust
//       region radius, so both must be positive with RHOEND no greater than
//       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
//       expected change to a variable, while RHOEND should indicate the
//       accuracy that is required in the final values of the variables. An
//       error return occurs if any of the differences XU[I]-XL[I], I=1,...,N,
//       is less than 2*RHOBEG.
//     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
//       amount of printing. Specifically, there is no output if IPRINT=0 and
//       there is output only at the return if IPRINT=1. Otherwise, each new
//       value of RHO is printed, with the best vector of variables so far and
//       the corresponding value of the objective function. Further, each new
//       value of F with its variables are output if IPRINT=3.
//     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
//     The array W will be used for working space. Its length must be at least
//       (NPT+5)*(NPT+N)+3*N*(N+5)/2.
//
//     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
//     F to the value of the objective function for the current values of the
//     variables X(1),X(2),...,X(N), which are generated automatically in a
//     way that satisfies the bounds given in XL and XU.
//
//     Return if the value of NPT is unacceptable.
//
      NP=N+1
      if (NPT < N+2 || NPT > ((N+2)*NP)/2) {
          PRINT 10
   10     FORMAT (/4X,'Return from BOBYQA because NPT is not in',
     1      ' the required interval')
          GO TO 40
      }
//
//     Partition the working space array, so that different parts of it can
//     be treated separately during the calculation of BOBYQB. The partition
//     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
//     space that is taken by the last array in the argument list of BOBYQB.
//
      NDIM=NPT+N
      IXB=1
      IXP=IXB+N
      IFV=IXP+N*NPT
      IXO=IFV+NPT
      IGO=IXO+N
      IHQ=IGO+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ISL=IZMAT+NPT*(NPT-NP)
      ISU=ISL+N
      IXN=ISU+N
      IXA=IXN+N
      ID=IXA+N
      IVL=ID+N
      IW=IVL+NDIM
//
//     Return if there is insufficient space between the bounds. Modify the
//     initial X if necessary in order to avoid conflicts between the bounds
//     and the construction of the first quadratic model. The lower and upper
//     bounds on moves from the updated X are set now, in the ISL and ISU
//     partitions of W, in order to provide useful and exact information about
//     components of X that become within distance RHOBEG from their bounds.
//
      ZERO=0.0D0
      DO 30 J=1,N
      TEMP=XU[J]-XL[J]
      if (TEMP < RHOBEG+RHOBEG) {
          PRINT 20
   20     FORMAT (/4X,'Return from BOBYQA because one of the',
     1      ' differences XU[I]-XL[I]'/6X,' is less than 2*RHOBEG.')
          GO TO 40
      }
      JSL=ISL+J-1
      JSU=JSL+N
      W(JSL)=XL[J]-X[J]
      W(JSU)=XU[J]-X[J]
      if (W(JSL) >= -RHOBEG) {
          if (W(JSL) >= ZERO) {
              X[J]=XL[J]
              W(JSL)=ZERO
              W(JSU)=TEMP
          ELSE
              X[J]=XL[J]+RHOBEG
              W(JSL)=-RHOBEG
              W(JSU)=Math.Max(XU[J]-X[J],RHOBEG)
          }
      ELSE if (W(JSU) <= RHOBEG) {
          if (W(JSU) <= ZERO) {
              X[J]=XU[J]
              W(JSL)=-TEMP
              W(JSU)=ZERO
          ELSE
              X[J]=XU[J]-RHOBEG
              W(JSL)=Math.Min(XL[J]-X[J],-RHOBEG)
              W(JSU)=RHOBEG
          }
      }
   30 CONTINUE
//
//     Make the call of BOBYQB.
//
      CALL BOBYQB (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),
     1  W(IXP),W(IFV),W(IXO),W(IGO),W(IHQ),W(IPQ),W(IBMAT),W(IZMAT),
     2  NDIM,W(ISL),W(ISU),W(IXN),W(IXA),W(ID),W(IVL),W(IW))
   40 RETURN
      END


      SUBROUTINE BOBYQB (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,
     1  MAXFUN,XBASE,XPT,FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,
     2  SL,SU,XNEW,XALT,D,VLAG,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),
     1  XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),
     2  SL(*),SU(*),XNEW(*),XALT(*),D(*),VLAG(*),W(*)
//
//     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN
//       are identical to the corresponding arguments in SUBROUTINE BOBYQA.
//     XBASE holds a shift of origin that should reduce the contributions
//       from rounding errors to values of the model and Lagrange functions.
//     XPT is a two-dimensional array that holds the coordinates of the
//       interpolation points relative to XBASE.
//     FVAL holds the values of F at the interpolation points.
//     XOPT is set to the displacement from XBASE of the trust region centre.
//     GOPT holds the gradient of the quadratic model at XBASE+XOPT.
//     HQ holds the explicit second derivatives of the quadratic model.
//     PQ contains the parameters of the implicit second derivatives of the
//       quadratic model.
//     BMAT holds the last N columns of H.
//     ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
//       this factorization being ZMAT times ZMAT^T, which provides both the
//       correct rank and positive semi-definiteness.
//     NDIM is the first dimension of BMAT and has the value NPT+N.
//     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
//       All the components of every XOPT are going to satisfy the bounds
//       SL[I] .LEQ. XOPT[I] .LEQ. SU[I], with appropriate equalities when
//       XOPT is on a constraint boundary.
//     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
//       vector of variables for the next call of CALFUN. XNEW also satisfies
//       the SL and SU constraints in the way that has just been mentioned.
//     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
//       in order to increase the denominator in the updating of UPDATE.
//     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
//     VLAG contains the values of the Lagrange functions at a new point X.
//       They are part of a product that requires VLAG to be of length NDIM.
//     W is a one-dimensional array that is used for working space. Its length
//       must be at least 3*NDIM = 3*(NPT+N).
//
//     Set some constants.
//
      HALF=0.5D0
      ONE=1.0D0
      TEN=10.0D0
      TENTH=0.1D0
      TWO=2.0D0
      ZERO=0.0D0
      NP=N+1
      NPTM=NPT-NP
      NH=(N*NP)/2
//
//     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
//     BMAT and ZMAT for the first iteration, with the corresponding values of
//     of NF and KOPT, which are the number of calls of CALFUN so far and the
//     index of the interpolation point at the trust region centre. Then the
//     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
//     less than NPT. GOPT will be updated if KOPT is different from KBASE.
//
      CALL PRELIM (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,XPT,
     1  FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT)
      XOPTSQ=ZERO
      DO 10 I=1,N
      XOPT[I]=XPT(KOPT,I)
   10 XOPTSQ=XOPTSQ+XOPT[I]**2
      FSAVE=FVAL(1)
      if (NF < NPT) {
          if (IPRINT > 0) PRINT 390
          goto 720
      }
      KBASE=1
//
//     Complete the settings that are required for the iterative procedure.
//
      RHO=RHOBEG
      DELTA=RHO
      NRESC=NF
      NTRITS=0
      DIFFA=ZERO
      DIFFB=ZERO
      ITEST=0
      NFSAV=NF
//
//     Update GOPT if necessary before the first iteration and after each
//     call of RESCUE that makes a call of CALFUN.
//
   20 if (KOPT != KBASE) {
          IH=0
          DO 30 J=1,N
          DO 30 I=1,J
          IH=IH+1
          if (I < J) GOPT[J]=GOPT[J]+HQ(IH)*XOPT[I]
   30     GOPT[I]=GOPT[I]+HQ(IH)*XOPT[J]
          if (NF > NPT) {
              DO 50 K=1,NPT
              TEMP=ZERO
              DO 40 J=1,N
   40         TEMP=TEMP+XPT(K,J)*XOPT[J]
              TEMP=PQ[K]*TEMP
              DO 50 I=1,N
   50         GOPT[I]=GOPT[I]+TEMP*XPT(K,I)
          }
      }
//
//     Generate the next point in the trust region that provides a small value
//     of the quadratic model subject to the constraints on the variables.
//     The integer NTRITS is set to the number "trust region" iterations that
//     have occurred since the last "alternative" iteration. If the length
//     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
//     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.
//
   60 CALL TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,XNEW,D,
     1  W,W(NP),W(NP+N),W(NP+2*N),W(NP+3*N),DSQ,CRVMIN)
      DNORM=Math.Min(DELTA,Math.Sqrt(DSQ))
      if (DNORM < HALF*RHO) {
          NTRITS=-1
          DISTSQ=(TEN*RHO)**2
          if (NF <= NFSAV+2) goto 650
//
//     The following choice between labels 650 and 680 depends on whether or
//     not our work with the current RHO seems to be complete. Either RHO is
//     decreased or termination occurs if the errors in the quadratic model at
//     the last three interpolation points compare favourably with predictions
//     of likely improvements to the model within distance HALF*RHO of XOPT.
//
          ERRBIG=Math.Max(DIFFA,DIFFB,DIFFC)
          FRHOSQ=0.125D0*RHO*RHO
          if (CRVMIN > ZERO && ERRBIG > FRHOSQ*CRVMIN)
     1       goto 650
          BDTOL=ERRBIG/RHO
          DO 80 J=1,N
          BDTEST=BDTOL
          if (XNEW[J] == SL[J]) BDTEST=W[J]
          if (XNEW[J] == SU[J]) BDTEST=-W[J]
          if (BDTEST < BDTOL) {
              CURV=HQ((J+J*J)/2)
              DO 70 K=1,NPT
   70         CURV=CURV+PQ[K]*XPT(K,J)**2
              BDTEST=BDTEST+HALF*CURV*RHO
              if (BDTEST < BDTOL) goto 650
          }
   80     CONTINUE
          goto 680
      }
      NTRITS=NTRITS+1
//
//     Severe cancellation is likely to occur if XOPT is too far from XBASE.
//     If the following test holds, then XBASE is shifted so that XOPT becomes
//     zero. The appropriate changes are made to BMAT and to the second
//     derivatives of the current model, beginning with the changes to BMAT
//     that do not depend on ZMAT. VLAG is used temporarily for working space.
//
   90 if (DSQ <= 1.0D-3*XOPTSQ) {
          FRACSQ=0.25D0*XOPTSQ
          SUMPQ=ZERO
          DO 110 K=1,NPT
          SUMPQ=SUMPQ+PQ[K]
          SUM=-HALF*XOPTSQ
          DO 100 I=1,N
  100     SUM=SUM+XPT(K,I)*XOPT[I]
          W(NPT+K)=SUM
          TEMP=FRACSQ-HALF*SUM
          DO 110 I=1,N
          W[I]=BMAT(K,I)
          VLAG[I]=SUM*XPT(K,I)+TEMP*XOPT[I]
          IP=NPT+I
          DO 110 J=1,I
  110     BMAT(IP,J)=BMAT(IP,J)+W[I]*VLAG[J]+VLAG[I]*W[J]
//
//     Then the revisions of BMAT that depend on ZMAT are calculated.
//
          DO 150 JJ=1,NPTM
          SUMZ=ZERO
          SUMW=ZERO
          DO 120 K=1,NPT
          SUMZ=SUMZ+ZMAT(K,JJ)
          VLAG[K]=W(NPT+K)*ZMAT(K,JJ)
  120     SUMW=SUMW+VLAG[K]
          DO 140 J=1,N
          SUM=(FRACSQ*SUMZ-HALF*SUMW)*XOPT[J]
          DO 130 K=1,NPT
  130     SUM=SUM+VLAG[K]*XPT(K,J)
          W[J]=SUM
          DO 140 K=1,NPT
  140     BMAT(K,J)=BMAT(K,J)+SUM*ZMAT(K,JJ)
          DO 150 I=1,N
          IP=I+NPT
          TEMP=W[I]
          DO 150 J=1,I
  150     BMAT(IP,J)=BMAT(IP,J)+TEMP*W[J]
//
//     The following instructions complete the shift, including the changes
//     to the second derivative parameters of the quadratic model.
//
          IH=0
          DO 170 J=1,N
          W[J]=-HALF*SUMPQ*XOPT[J]
          DO 160 K=1,NPT
          W[J]=W[J]+PQ[K]*XPT(K,J)
  160     XPT(K,J)=XPT(K,J)-XOPT[J]
          DO 170 I=1,J
          IH=IH+1
          HQ(IH)=HQ(IH)+W[I]*XOPT[J]+XOPT[I]*W[J]
  170     BMAT(NPT+I,J)=BMAT(NPT+J,I)
          DO 180 I=1,N
          XBASE[I]=XBASE[I]+XOPT[I]
          XNEW[I]=XNEW[I]-XOPT[I]
          SL[I]=SL[I]-XOPT[I]
          SU[I]=SU[I]-XOPT[I]
  180     XOPT[I]=ZERO
          XOPTSQ=ZERO
      }
      if (NTRITS == 0) goto 210
      goto 230
//
//     XBASE is also moved to XOPT by a call of RESCUE. This calculation is
//     more expensive than the previous shift, because new matrices BMAT and
//     ZMAT are generated from scratch, which may include the replacement of
//     interpolation points whose positions seem to be causing near linear
//     dependence in the interpolation conditions. Therefore RESCUE is called
//     only if rounding errors have reduced by at least a factor of two the
//     denominator of the formula for updating the H matrix. It provides a
//     useful safeguard, but is not invoked in most applications of BOBYQA.
//
  190 NFSAV=NF
      KBASE=KOPT
      CALL RESCUE (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,FVAL,
     1  XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,KOPT,
     2  VLAG,W,W(N+NP),W(NDIM+NP))
//
//     XOPT is updated now in case the branch below to label 720 is taken.
//     Any updating of GOPT occurs after the branch below to label 20, which
//     leads to a trust region iteration as does the branch to label 60.
//
      XOPTSQ=ZERO
      if (KOPT != KBASE) {
          DO 200 I=1,N
          XOPT[I]=XPT(KOPT,I)
  200     XOPTSQ=XOPTSQ+XOPT[I]**2
      }
      if (NF < 0) {
          NF=MAXFUN
          if (IPRINT > 0) PRINT 390
          goto 720
      }
      NRESC=NF
      if (NFSAV < NF) {
          NFSAV=NF
          goto 20
      }
      if (NTRITS > 0) goto 60
//
//     Pick two alternative vectors of variables, relative to XBASE, that
//     are suitable as new positions of the KNEW-th interpolation point.
//     Firstly, XNEW is set to the point on a line through XOPT and another
//     interpolation point that minimizes the predicted value of the next
//     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
//     and SU bounds. Secondly, XALT is set to the best feasible point on
//     a constrained version of the Cauchy step of the KNEW-th Lagrange
//     function, the corresponding value of the square of this function
//     being returned in CAUCHY. The choice between these alternatives is
//     going to be made when the denominator is calculated.
//
  210 CALL ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,
     1  KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,W,W(NP),W(NDIM+1))
      DO 220 I=1,N
  220 D[I]=XNEW[I]-XOPT[I]
//
//     Calculate VLAG and BETA for the current choice of D. The scalar
//     product of D with XPT(K,.) is going to be held in W(NPT+K) for
//     use when VQUAD is calculated.
//
  230 DO 250 K=1,NPT
      SUMA=ZERO
      SUMB=ZERO
      SUM=ZERO
      DO 240 J=1,N
      SUMA=SUMA+XPT(K,J)*D[J]
      SUMB=SUMB+XPT(K,J)*XOPT[J]
  240 SUM=SUM+BMAT(K,J)*D[J]
      W[K]=SUMA*(HALF*SUMA+SUMB)
      VLAG[K]=SUM
  250 W(NPT+K)=SUMA
      BETA=ZERO
      DO 270 JJ=1,NPTM
      SUM=ZERO
      DO 260 K=1,NPT
  260 SUM=SUM+ZMAT(K,JJ)*W[K]
      BETA=BETA-SUM*SUM
      DO 270 K=1,NPT
  270 VLAG[K]=VLAG[K]+SUM*ZMAT(K,JJ)
      DSQ=ZERO
      BSUM=ZERO
      DX=ZERO
      DO 300 J=1,N
      DSQ=DSQ+D[J]**2
      SUM=ZERO
      DO 280 K=1,NPT
  280 SUM=SUM+W[K]*BMAT(K,J)
      BSUM=BSUM+SUM*D[J]
      JP=NPT+J
      DO 290 I=1,N
  290 SUM=SUM+BMAT(JP,I)*D[I]
      VLAG(JP)=SUM
      BSUM=BSUM+SUM*D[J]
  300 DX=DX+D[J]*XOPT[J]
      BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
//
//     If NTRITS is zero, the denominator may be increased by replacing
//     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
//     rounding errors have damaged the chosen denominator.
//
      if (NTRITS == 0) {
          DENOM=VLAG(KNEW)**2+ALPHA*BETA
          if (DENOM < CAUCHY && CAUCHY > ZERO) {
              DO 310 I=1,N
              XNEW[I]=XALT[I]
  310         D[I]=XNEW[I]-XOPT[I]
              CAUCHY=ZERO
              GO TO 230
          }
          if (DENOM <= HALF*VLAG(KNEW)**2) {
              if (NF > NRESC) goto 190
              if (IPRINT > 0) PRINT 320
  320         FORMAT (/5X,'Return from BOBYQA because of much',
     1          ' cancellation in a denominator.')
              goto 720
          }
//
//     Alternatively, if NTRITS is positive, then set KNEW to the index of
//     the next interpolation point to be deleted to make room for a trust
//     region step. Again RESCUE may be called if rounding errors have damaged
//     the chosen denominator, which is the reason for attempting to select
//     KNEW before calculating the next value of the objective function.
//
      ELSE
          DELSQ=DELTA*DELTA
          SCADEN=ZERO
          BIGLSQ=ZERO
          KNEW=0
          DO 350 K=1,NPT
          if (K == KOPT) goto 350
          HDIAG=ZERO
          DO 330 JJ=1,NPTM
  330     HDIAG=HDIAG+ZMAT(K,JJ)**2
          DEN=BETA*HDIAG+VLAG[K]**2
          DISTSQ=ZERO
          DO 340 J=1,N
  340     DISTSQ=DISTSQ+(XPT(K,J)-XOPT[J])**2
          TEMP=Math.Max(ONE,(DISTSQ/DELSQ)**2)
          if (TEMP*DEN > SCADEN) {
              SCADEN=TEMP*DEN
              KNEW=K
              DENOM=DEN
          }
          BIGLSQ=Math.Max(BIGLSQ,TEMP*VLAG[K]**2)
  350     CONTINUE
          if (SCADEN <= HALF*BIGLSQ) {
              if (NF > NRESC) goto 190
              if (IPRINT > 0) PRINT 320
              goto 720
          }
      }
//
//     Put the variables for the next calculation of the objective function
//       in XNEW, with any adjustments for the bounds.
//
//
//     Calculate the value of the objective function at XBASE+XNEW, unless
//       the limit on the number of calculations of F has been reached.
//
  360 DO 380 I=1,N
      X[I]=Math.Min(Math.Max(XL[I],XBASE[I]+XNEW[I]),XU[I])
      if (XNEW[I] == SL[I]) X[I]=XL[I]
      if (XNEW[I] == SU[I]) X[I]=XU[I]
  380 CONTINUE
      if (NF >= MAXFUN) {
          if (IPRINT > 0) PRINT 390
  390     FORMAT (/4X,'Return from BOBYQA because CALFUN has been',
     1      ' called MAXFUN times.')
          goto 720
      }
      NF=NF+1
      CALL CALFUN (N,X,F)
      if (IPRINT == 3) {
          PRINT 400, NF,F,(X[I],I=1,N)
  400      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      }
      if (NTRITS == -1) {
          FSAVE=F
          goto 720
      }
//
//     Use the quadratic model to predict the change in F due to the step D,
//       and set DIFF to the error of this prediction.
//
      FOPT=FVAL(KOPT)
      VQUAD=ZERO
      IH=0
      DO 410 J=1,N
      VQUAD=VQUAD+D[J]*GOPT[J]
      DO 410 I=1,J
      IH=IH+1
      TEMP=D[I]*D[J]
      if (I == J) TEMP=HALF*TEMP
  410 VQUAD=VQUAD+HQ(IH)*TEMP
      DO 420 K=1,NPT
  420 VQUAD=VQUAD+HALF*PQ[K]*W(NPT+K)**2
      DIFF=F-FOPT-VQUAD
      DIFFC=DIFFB
      DIFFB=DIFFA
      DIFFA=DABS(DIFF)
      if (DNORM > RHO) NFSAV=NF
//
//     Pick the next value of DELTA after a trust region step.
//
      if (NTRITS > 0) {
          if (VQUAD >= ZERO) {
              if (IPRINT > 0) PRINT 430
  430         FORMAT (/4X,'Return from BOBYQA because a trust',
     1          ' region step has failed to reduce Q.')
              goto 720
          }
          RATIO=(F-FOPT)/VQUAD
          if (RATIO <= TENTH) {
              DELTA=Math.Min(HALF*DELTA,DNORM)
          ELSE if (RATIO. LE. 0.7D0) {
              DELTA=Math.Max(HALF*DELTA,DNORM)
          ELSE
              DELTA=Math.Max(HALF*DELTA,DNORM+DNORM)
          }
          if (DELTA <= 1.5D0*RHO) DELTA=RHO
//
//     Recalculate KNEW and DENOM if the new F is less than FOPT.
//
          if (F < FOPT) {
              KSAV=KNEW
              DENSAV=DENOM
              DELSQ=DELTA*DELTA
              SCADEN=ZERO
              BIGLSQ=ZERO
              KNEW=0
              DO 460 K=1,NPT
              HDIAG=ZERO
              DO 440 JJ=1,NPTM
  440         HDIAG=HDIAG+ZMAT(K,JJ)**2
              DEN=BETA*HDIAG+VLAG[K]**2
              DISTSQ=ZERO
              DO 450 J=1,N
  450         DISTSQ=DISTSQ+(XPT(K,J)-XNEW[J])**2
              TEMP=Math.Max(ONE,(DISTSQ/DELSQ)**2)
              if (TEMP*DEN > SCADEN) {
                  SCADEN=TEMP*DEN
                  KNEW=K
                  DENOM=DEN
              }
  460         BIGLSQ=Math.Max(BIGLSQ,TEMP*VLAG[K]**2)
              if (SCADEN <= HALF*BIGLSQ) {
                  KNEW=KSAV
                  DENOM=DENSAV
              }
          }
      }
//
//     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
//     moved. Also update the second derivative terms of the model.
//
      CALL UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
      IH=0
      PQOLD=PQ(KNEW)
      PQ(KNEW)=ZERO
      DO 470 I=1,N
      TEMP=PQOLD*XPT(KNEW,I)
      DO 470 J=1,I
      IH=IH+1
  470 HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
      DO 480 JJ=1,NPTM
      TEMP=DIFF*ZMAT(KNEW,JJ)
      DO 480 K=1,NPT
  480 PQ[K]=PQ[K]+TEMP*ZMAT(K,JJ)
//
//     Include the new interpolation point, and make the changes to GOPT at
//     the old XOPT that are caused by the updating of the quadratic model.
//
      FVAL(KNEW)=F
      DO 490 I=1,N
      XPT(KNEW,I)=XNEW[I]
  490 W[I]=BMAT(KNEW,I)
      DO 520 K=1,NPT
      SUMA=ZERO
      DO 500 JJ=1,NPTM
  500 SUMA=SUMA+ZMAT(KNEW,JJ)*ZMAT(K,JJ)
      SUMB=ZERO
      DO 510 J=1,N
  510 SUMB=SUMB+XPT(K,J)*XOPT[J]
      TEMP=SUMA*SUMB
      DO 520 I=1,N
  520 W[I]=W[I]+TEMP*XPT(K,I)
      DO 530 I=1,N
  530 GOPT[I]=GOPT[I]+DIFF*W[I]
//
//     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
//
      if (F < FOPT) {
          KOPT=KNEW
          XOPTSQ=ZERO
          IH=0
          DO 540 J=1,N
          XOPT[J]=XNEW[J]
          XOPTSQ=XOPTSQ+XOPT[J]**2
          DO 540 I=1,J
          IH=IH+1
          if (I < J) GOPT[J]=GOPT[J]+HQ(IH)*D[I]
  540     GOPT[I]=GOPT[I]+HQ(IH)*D[J]
          DO 560 K=1,NPT
          TEMP=ZERO
          DO 550 J=1,N
  550     TEMP=TEMP+XPT(K,J)*D[J]
          TEMP=PQ[K]*TEMP
          DO 560 I=1,N
  560     GOPT[I]=GOPT[I]+TEMP*XPT(K,I)
      }
//
//     Calculate the parameters of the least Frobenius norm interpolant to
//     the current data, the gradient of this interpolant at XOPT being put
//     into VLAG(NPT+I), I=1,2,...,N.
//
      if (NTRITS > 0) {
          DO 570 K=1,NPT
          VLAG[K]=FVAL[K]-FVAL(KOPT)
  570     W[K]=ZERO
          DO 590 J=1,NPTM
          SUM=ZERO
          DO 580 K=1,NPT
  580     SUM=SUM+ZMAT(K,J)*VLAG[K]
          DO 590 K=1,NPT
  590     W[K]=W[K]+SUM*ZMAT(K,J)
          DO 610 K=1,NPT
          SUM=ZERO
          DO 600 J=1,N
  600     SUM=SUM+XPT(K,J)*XOPT[J]
          W(K+NPT)=W[K]
  610     W[K]=SUM*W[K]
          GQSQ=ZERO
          GISQ=ZERO
          DO 630 I=1,N
          SUM=ZERO
          DO 620 K=1,NPT
  620     SUM=SUM+BMAT(K,I)*VLAG[K]+XPT(K,I)*W[K]
          if (XOPT[I] == SL[I]) {
              GQSQ=GQSQ+Math.Min(ZERO,GOPT[I])**2
              GISQ=GISQ+Math.Min(ZERO,SUM)**2
          ELSE if (XOPT[I] == SU[I]) {
              GQSQ=GQSQ+Math.Max(ZERO,GOPT[I])**2
              GISQ=GISQ+Math.Max(ZERO,SUM)**2
          ELSE
              GQSQ=GQSQ+GOPT[I]**2
              GISQ=GISQ+SUM*SUM
          }
  630     VLAG(NPT+I)=SUM
//
//     Test whether to replace the new quadratic model by the least Frobenius
//     norm interpolant, making the replacement if the test is satisfied.
//
          ITEST=ITEST+1
          if (GQSQ < TEN*GISQ) ITEST=0
          if (ITEST >= 3) {
              DO 640 I=1,MAX0(NPT,NH)
              if (I <= N) GOPT[I]=VLAG(NPT+I)
              if (I <= NPT) PQ[I]=W(NPT+I)
              if (I <= NH) HQ[I]=ZERO
              ITEST=0
  640         CONTINUE
          }
      }
//
//     If a trust region step has provided a sufficient decrease in F, then
//     branch for another trust region calculation. The case NTRITS=0 occurs
//     when the new interpolation point was reached by an alternative step.
//
      if (NTRITS == 0) goto 60
      if (F <= FOPT+TENTH*VQUAD) goto 60
//
//     Alternatively, find out if the interpolation points are close enough
//       to the best point so far.
//
      DISTSQ=Math.Max((TWO*DELTA)**2,(TEN*RHO)**2)
  650 KNEW=0
      DO 670 K=1,NPT
      SUM=ZERO
      DO 660 J=1,N
  660 SUM=SUM+(XPT(K,J)-XOPT[J])**2
      if (SUM > DISTSQ) {
          KNEW=K
          DISTSQ=SUM
      }
  670 CONTINUE
//
//     If KNEW is positive, then ALTMOV finds alternative new positions for
//     the KNEW-th interpolation point within distance ADELT of XOPT. It is
//     reached via label 90. Otherwise, there is a branch to label 60 for
//     another trust region iteration, unless the calculations with the
//     current RHO are complete.
//
      if (KNEW > 0) {
          DIST=Math.Sqrt(DISTSQ)
          if (NTRITS == -1) {
              DELTA=Math.Min(TENTH*DELTA,HALF*DIST)
              if (DELTA <= 1.5D0*RHO) DELTA=RHO
          }
          NTRITS=0
          ADELT=Math.Max(Math.Min(TENTH*DIST,DELTA),RHO)
          DSQ=ADELT*ADELT
          goto 90
      }
      if (NTRITS == -1) goto 680
      if (RATIO > ZERO) goto 60
      if (Math.Max(DELTA,DNORM) > RHO) goto 60
//
//     The calculations with the current value of RHO are complete. Pick the
//       next values of RHO and DELTA.
//
  680 if (RHO > RHOEND) {
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          if (RATIO <= 16.0D0) {
              RHO=RHOEND
          ELSE if (RATIO <= 250.0D0) {
              RHO=Math.Sqrt(RATIO)*RHOEND
          ELSE
              RHO=TENTH*RHO
          }
          DELTA=Math.Max(DELTA,RHO)
          if (IPRINT >= 2) {
              if (IPRINT >= 3) PRINT 690
  690         FORMAT (5X)
              PRINT 700, RHO,NF
  700         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
     1          ' function values =',I6)
              PRINT 710, FVAL(KOPT),(XBASE[I]+XOPT[I],I=1,N)
  710         FORMAT (4X,'Least value of F =',1PD23.15,9X,
     1          'The corresponding X is:'/(2X,5D15.6))
          }
          NTRITS=0
          NFSAV=NF
          goto 60
      }
//
//     Return from the calculation, after another Newton-Raphson step, if
//       it is too short to have been tried before.
//
      if (NTRITS == -1) goto 360
  720 if (FVAL(KOPT) <= FSAVE) {
          DO 730 I=1,N
          X[I]=Math.Min(Math.Max(XL[I],XBASE[I]+XOPT[I]),XU[I])
          if (XOPT[I] == SL[I]) X[I]=XL[I]
          if (XOPT[I] == SU[I]) X[I]=XU[I]
  730     CONTINUE
          F=FVAL(KOPT)
      }
      if (IPRINT >= 1) {
          PRINT 740, NF
  740     FORMAT (/4X,'At the return from BOBYQA',5X,
     1      'Number of function values =',I6)
          PRINT 710, F,(X[I],I=1,N)
      }
      RETURN
      END


      SUBROUTINE ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,
     1  KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,GLAG,HCOL,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XPT(NPT,*),XOPT(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),
     1  SU(*),XNEW(*),XALT(*),GLAG(*),HCOL(*),W(*)
//
//     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
//       the same meanings as the corresponding arguments of BOBYQB.
//     KOPT is the index of the optimal interpolation point.
//     KNEW is the index of the interpolation point that is going to be moved.
//     ADELT is the current trust region bound.
//     XNEW will be set to a suitable new position for the interpolation point
//       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
//       bounds and it should provide a large denominator in the next call of
//       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
//       straight lines through XOPT and another interpolation point.
//     XALT also provides a large value of the modulus of the KNEW-th Lagrange
//       function subject to the constraints that have been mentioned, its main
//       difference from XNEW being that XALT-XOPT is a constrained version of
//       the Cauchy step within the trust region. An exception is that XALT is
//       not calculated if all components of GLAG (see below) are zero.
//     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
//     CAUCHY will be set to the square of the KNEW-th Lagrange function at
//       the step XALT-XOPT from XOPT for the vector XALT that is returned,
//       except that CAUCHY is set to zero if XALT is not calculated.
//     GLAG is a working space vector of length N for the gradient of the
//       KNEW-th Lagrange function at XOPT.
//     HCOL is a working space vector of length NPT for the second derivative
//       coefficients of the KNEW-th Lagrange function.
//     W is a working space vector of length 2N that is going to hold the
//       constrained Cauchy step from XOPT of the Lagrange function, followed
//       by the downhill version of XALT when the uphill step is calculated.
//
//     Set the first NPT components of W to the leading elements of the
//     KNEW-th column of the H matrix.
//
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      CONST=ONE+Math.Sqrt(2.0D0)
      DO 10 K=1,NPT
   10 HCOL[K]=ZERO
      DO 20 J=1,NPT-N-1
      TEMP=ZMAT(KNEW,J)
      DO 20 K=1,NPT
   20 HCOL[K]=HCOL[K]+TEMP*ZMAT(K,J)
      ALPHA=HCOL(KNEW)
      HA=HALF*ALPHA
//
//     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
//
      DO 30 I=1,N
   30 GLAG[I]=BMAT(KNEW,I)
      DO 50 K=1,NPT
      TEMP=ZERO
      DO 40 J=1,N
   40 TEMP=TEMP+XPT(K,J)*XOPT[J]
      TEMP=HCOL[K]*TEMP
      DO 50 I=1,N
   50 GLAG[I]=GLAG[I]+TEMP*XPT(K,I)
//
//     Search for a large denominator along the straight lines through XOPT
//     and another interpolation point. SLBD and SUBD will be lower and upper
//     bounds on the step along each of these lines in turn. PREDSQ will be
//     set to the square of the predicted denominator for each line. PRESAV
//     will be set to the largest admissible value of PREDSQ that occurs.
//
      PRESAV=ZERO
      DO 80 K=1,NPT
      if (K == KOPT) goto 80
      DDERIV=ZERO
      DISTSQ=ZERO
      DO 60 I=1,N
      TEMP=XPT(K,I)-XOPT[I]
      DDERIV=DDERIV+GLAG[I]*TEMP
   60 DISTSQ=DISTSQ+TEMP*TEMP
      SUBD=ADELT/Math.Sqrt(DISTSQ)
      SLBD=-SUBD
      ILBD=0
      IUBD=0
      SUMIN=Math.Min(ONE,SUBD)
//
//     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
//
      DO 70 I=1,N
      TEMP=XPT(K,I)-XOPT[I]
      if (TEMP > ZERO) {
          if (SLBD*TEMP < SL[I]-XOPT[I]) {
              SLBD=(SL[I]-XOPT[I])/TEMP
              ILBD=-I
          }
          if (SUBD*TEMP > SU[I]-XOPT[I]) {
              SUBD=Math.Max(SUMIN,(SU[I]-XOPT[I])/TEMP)
              IUBD=I
          }
      ELSE if (TEMP < ZERO) {
          if (SLBD*TEMP > SU[I]-XOPT[I]) {
              SLBD=(SU[I]-XOPT[I])/TEMP
              ILBD=I
          }
          if (SUBD*TEMP < SL[I]-XOPT[I]) {
              SUBD=Math.Max(SUMIN,(SL[I]-XOPT[I])/TEMP)
              IUBD=-I
          }
      }
   70 CONTINUE
//
//     Seek a large modulus of the KNEW-th Lagrange function when the index
//     of the other interpolation point on the line through XOPT is KNEW.
//
      if (K == KNEW) {
          DIFF=DDERIV-ONE
          STEP=SLBD
          VLAG=SLBD*(DDERIV-SLBD*DIFF)
          ISBD=ILBD
          TEMP=SUBD*(DDERIV-SUBD*DIFF)
          if (DABS(TEMP) > DABS(VLAG)) {
              STEP=SUBD
              VLAG=TEMP
              ISBD=IUBD
          }
          TEMPD=HALF*DDERIV
          TEMPA=TEMPD-DIFF*SLBD
          TEMPB=TEMPD-DIFF*SUBD
          if (TEMPA*TEMPB < ZERO) {
              TEMP=TEMPD*TEMPD/DIFF
              if (DABS(TEMP) > DABS(VLAG)) {
                  STEP=TEMPD/DIFF
                  VLAG=TEMP
                  ISBD=0
              }
          }
//
//     Search along each of the other lines through XOPT and another point.
//
      ELSE
          STEP=SLBD
          VLAG=SLBD*(ONE-SLBD)
          ISBD=ILBD
          TEMP=SUBD*(ONE-SUBD)
          if (DABS(TEMP) > DABS(VLAG)) {
              STEP=SUBD
              VLAG=TEMP
              ISBD=IUBD
          }
          if (SUBD > HALF) {
              if (DABS(VLAG) < 0.25D0) {
                  STEP=HALF
                  VLAG=0.25D0
                  ISBD=0
              }
          }
          VLAG=VLAG*DDERIV
      }
//
//     Calculate PREDSQ for the current line search and maintain PRESAV.
//
      TEMP=STEP*(ONE-STEP)*DISTSQ
      PREDSQ=VLAG*VLAG*(VLAG*VLAG+HA*TEMP*TEMP)
      if (PREDSQ > PRESAV) {
          PRESAV=PREDSQ
          KSAV=K
          STPSAV=STEP
          IBDSAV=ISBD
      }
   80 CONTINUE
//
//     Construct XNEW in a way that satisfies the bound constraints exactly.
//
      DO 90 I=1,N
      TEMP=XOPT[I]+STPSAV*(XPT(KSAV,I)-XOPT[I])
   90 XNEW[I]=Math.Max(SL[I],Math.Min(SU[I],TEMP))
      if (IBDSAV < 0) XNEW(-IBDSAV)=SL(-IBDSAV)
      if (IBDSAV > 0) XNEW(IBDSAV)=SU(IBDSAV)
//
//     Prepare for the iterative method that assembles the constrained Cauchy
//     step in W. The sum of squares of the fixed components of W is formed in
//     WFIXSQ, and the free components of W are set to BIGSTP.
//
      BIGSTP=ADELT+ADELT
      IFLAG=0
  100 WFIXSQ=ZERO
      GGFREE=ZERO
      DO 110 I=1,N
      W[I]=ZERO
      TEMPA=Math.Min(XOPT[I]-SL[I],GLAG[I])
      TEMPB=Math.Max(XOPT[I]-SU[I],GLAG[I])
      if (TEMPA > ZERO || TEMPB < ZERO) {
          W[I]=BIGSTP
          GGFREE=GGFREE+GLAG[I]**2
      }
  110 CONTINUE
      if (GGFREE == ZERO) {
          CAUCHY=ZERO
          goto 200
      }
//
//     Investigate whether more components of W can be fixed.
//
  120 TEMP=ADELT*ADELT-WFIXSQ
      if (TEMP > ZERO) {
          WSQSAV=WFIXSQ
          STEP=Math.Sqrt(TEMP/GGFREE)
          GGFREE=ZERO
          DO 130 I=1,N
          if (W[I] == BIGSTP) {
              TEMP=XOPT[I]-STEP*GLAG[I]
              if (TEMP <= SL[I]) {
                  W[I]=SL[I]-XOPT[I]
                  WFIXSQ=WFIXSQ+W[I]**2
              ELSE if (TEMP >= SU[I]) {
                  W[I]=SU[I]-XOPT[I]
                  WFIXSQ=WFIXSQ+W[I]**2
              ELSE
                  GGFREE=GGFREE+GLAG[I]**2
              }
          }
  130     CONTINUE
          if (WFIXSQ > WSQSAV && GGFREE > ZERO) goto 120
      }
//
//     Set the remaining free components of W and all components of XALT,
//     except that W may be scaled later.
//
      GW=ZERO
      DO 140 I=1,N
      if (W[I] == BIGSTP) {
          W[I]=-STEP*GLAG[I]
          XALT[I]=Math.Max(SL[I],Math.Min(SU[I],XOPT[I]+W[I]))
      ELSE if (W[I] == ZERO) {
          XALT[I]=XOPT[I]
      ELSE if (GLAG[I] > ZERO) {
          XALT[I]=SL[I]
      ELSE
          XALT[I]=SU[I]
      }
  140 GW=GW+GLAG[I]*W[I]
//
//     Set CURV to the curvature of the KNEW-th Lagrange function along W.
//     Scale W by a factor less than one if that can reduce the modulus of
//     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
//     the square of this function.
//
      CURV=ZERO
      DO 160 K=1,NPT
      TEMP=ZERO
      DO 150 J=1,N
  150 TEMP=TEMP+XPT(K,J)*W[J]
  160 CURV=CURV+HCOL[K]*TEMP*TEMP
      if (IFLAG == 1) CURV=-CURV
      if (CURV > -GW && CURV < -CONST*GW) {
          SCALE=-GW/CURV
          DO 170 I=1,N
          TEMP=XOPT[I]+SCALE*W[I]
  170     XALT[I]=Math.Max(SL[I],Math.Min(SU[I],TEMP))
          CAUCHY=(HALF*GW*SCALE)**2
      ELSE
          CAUCHY=(GW+HALF*CURV)**2
      }
//
//     If IFLAG is zero, then XALT is calculated as before after reversing
//     the sign of GLAG. Thus two XALT vectors become available. The one that
//     is chosen is the one that gives the larger value of CAUCHY.
//
      if (IFLAG == 0) {
          DO 180 I=1,N
          GLAG[I]=-GLAG[I]
  180     W(N+I)=XALT[I]
          CSAVE=CAUCHY
          IFLAG=1
          goto 100
      }
      if (CSAVE > CAUCHY) {
          DO 190 I=1,N
  190     XALT[I]=W(N+I)
          CAUCHY=CSAVE
      }
  200 RETURN
      END


      SUBROUTINE PRELIM (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,
     1  XPT,FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),GOPT(*),
     1  HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*)
//
//     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
//       same as the corresponding arguments in SUBROUTINE BOBYQA.
//     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
//       are the same as the corresponding arguments in BOBYQB, the elements
//       of SL and SU being set in BOBYQA.
//     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
//       it is set by PRELIM to the gradient of the quadratic model at XBASE.
//       If XOPT is nonzero, BOBYQB will change it to its usual value later.
//     NF is maintaned as the number of calls of CALFUN so far.
//     KOPT will be such that the least calculated value of F so far is at
//       the point XPT(KOPT,.)+XBASE in the space of the variables.
//
//     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
//     BMAT and ZMAT for the first iteration, and it maintains the values of
//     NF and KOPT. The vector X is also changed by PRELIM.
//
//     Set some constants.
//
      HALF=0.5D0
      ONE=1.0D0
      TWO=2.0D0
      ZERO=0.0D0
      RHOSQ=RHOBEG*RHOBEG
      RECIP=ONE/RHOSQ
      NP=N+1
//
//     Set XBASE to the initial vector of variables, and set the initial
//     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
//
      DO 20 J=1,N
      XBASE[J]=X[J]
      DO 10 K=1,NPT
   10 XPT(K,J)=ZERO
      DO 20 I=1,NDIM
   20 BMAT(I,J)=ZERO
      DO 30 IH=1,(N*NP)/2
   30 HQ(IH)=ZERO
      DO 40 K=1,NPT
      PQ[K]=ZERO
      DO 40 J=1,NPT-NP
   40 ZMAT(K,J)=ZERO
//
//     Begin the initialization procedure. NF becomes one more than the number
//     of function values so far. The coordinates of the displacement of the
//     next initial interpolation point from XBASE are set in XPT(NF+1,.).
//
      NF=0
   50 NFM=NF
      NFX=NF-N
      NF=NF+1
      if (NFM <= 2*N) {
          if (NFM >= 1 && NFM <= N) {
              STEPA=RHOBEG
              if (SU(NFM) == ZERO) STEPA=-STEPA
              XPT(NF,NFM)=STEPA
          ELSE if (NFM > N) {
              STEPA=XPT(NF-N,NFX)
              STEPB=-RHOBEG
              if (SL(NFX) == ZERO) STEPB=Math.Min(TWO*RHOBEG,SU(NFX))
              if (SU(NFX) == ZERO) STEPB=Math.Max(-TWO*RHOBEG,SL(NFX))
              XPT(NF,NFX)=STEPB
          }
      ELSE
          ITEMP=(NFM-NP)/N
          JPT=NFM-ITEMP*N-N
          IPT=JPT+ITEMP
          if (IPT > N) {
              ITEMP=JPT
              JPT=IPT-N
              IPT=ITEMP
          }
          XPT(NF,IPT)=XPT(IPT+1,IPT)
          XPT(NF,JPT)=XPT(JPT+1,JPT)
      }
//
//     Calculate the next value of F. The least function value so far and
//     its index are required.
//
      DO 60 J=1,N
      X[J]=Math.Min(Math.Max(XL[J],XBASE[J]+XPT(NF,J)),XU[J])
      if (XPT(NF,J) == SL[J]) X[J]=XL[J]
      if (XPT(NF,J) == SU[J]) X[J]=XU[J]
   60 CONTINUE
      CALL CALFUN (N,X,F)
      if (IPRINT == 3) {
          PRINT 70, NF,F,(X[I],I=1,N)
   70      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      }
      FVAL(NF)=F
      if (NF == 1) {
          FBEG=F
          KOPT=1
      ELSE if (F < FVAL(KOPT)) {
          KOPT=NF
      }
//
//     Set the nonzero initial elements of BMAT and the quadratic model in the
//     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
//     of the NF-th and (NF-N)-th interpolation points may be switched, in
//     order that the function value at the first of them contributes to the
//     off-diagonal second derivative terms of the initial quadratic model.
//
      if (NF <= 2*N+1) {
          if (NF >= 2 && NF <= N+1) {
              GOPT(NFM)=(F-FBEG)/STEPA
              if (NPT < NF+N) {
                  BMAT(1,NFM)=-ONE/STEPA
                  BMAT(NF,NFM)=ONE/STEPA
                  BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
              }
          ELSE if (NF >= N+2) {
              IH=(NFX*(NFX+1))/2
              TEMP=(F-FBEG)/STEPB
              DIFF=STEPB-STEPA
              HQ(IH)=TWO*(TEMP-GOPT(NFX))/DIFF
              GOPT(NFX)=(GOPT(NFX)*STEPB-TEMP*STEPA)/DIFF
              if (STEPA*STEPB < ZERO) {
                  if (F < FVAL(NF-N)) {
                      FVAL(NF)=FVAL(NF-N)
                      FVAL(NF-N)=F
                      if (KOPT == NF) KOPT=NF-N
                      XPT(NF-N,NFX)=STEPB
                      XPT(NF,NFX)=STEPA
                  }
              }
              BMAT(1,NFX)=-(STEPA+STEPB)/(STEPA*STEPB)
              BMAT(NF,NFX)=-HALF/XPT(NF-N,NFX)
              BMAT(NF-N,NFX)=-BMAT(1,NFX)-BMAT(NF,NFX)
              ZMAT(1,NFX)=Math.Sqrt(TWO)/(STEPA*STEPB)
              ZMAT(NF,NFX)=Math.Sqrt(HALF)/RHOSQ
              ZMAT(NF-N,NFX)=-ZMAT(1,NFX)-ZMAT(NF,NFX)
          }
//
//     Set the off-diagonal second derivatives of the Lagrange functions and
//     the initial quadratic model.
//
      ELSE
          IH=(IPT*(IPT-1))/2+JPT
          ZMAT(1,NFX)=RECIP
          ZMAT(NF,NFX)=RECIP
          ZMAT(IPT+1,NFX)=-RECIP
          ZMAT(JPT+1,NFX)=-RECIP
          TEMP=XPT(NF,IPT)*XPT(NF,JPT)
          HQ(IH)=(FBEG-FVAL(IPT+1)-FVAL(JPT+1)+F)/TEMP
      }
      if (NF < NPT && NF < MAXFUN) goto 50
      RETURN
      END


      SUBROUTINE RESCUE (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,
     1  FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,
     2  KOPT,VLAG,PTSAUX,PTSID,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),XOPT(*),
     1  GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*),
     2  VLAG(*),PTSAUX(2,*),PTSID(*),W(*)
//
//     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
//       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
//       the corresponding arguments of BOBYQB on the entry to RESCUE.
//     NF is maintained as the number of calls of CALFUN so far, except that
//       NF is set to -1 if the value of MAXFUN prevents further progress.
//     KOPT is maintained so that FVAL(KOPT) is the least calculated function
//       value. Its correct value must be given on entry. It is updated if a
//       new least function value is found, but the corresponding changes to
//       XOPT and GOPT have to be made later by the calling program.
//     DELTA is the current trust region radius.
//     VLAG is a working space vector that will be used for the values of the
//       provisional Lagrange functions at each of the interpolation points.
//       They are part of a product that requires VLAG to be of length NDIM.
//     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and
//       PTSAUX(2,J) specify the two positions of provisional interpolation
//       points when a nonzero step is taken along e_J (the J-th coordinate
//       direction) through XBASE+XOPT, as specified below. Usually these
//       steps have length DELTA, but other lengths are chosen if necessary
//       in order to satisfy the given bounds on the variables.
//     PTSID is also a working space array. It has NPT components that denote
//       provisional new positions of the original interpolation points, in
//       case changes are needed to restore the linear independence of the
//       interpolation conditions. The K-th point is a candidate for change
//       if and only if PTSID[K] is nonzero. In this case let p and q be the
//       integer parts of PTSID[K] and (PTSID[K]-p) multiplied by N+1. If p
//       and q are both positive, the step from XBASE+XOPT to the new K-th
//       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
//       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
//       p=0, respectively.
//     The first NDIM+NPT elements of the array W are used for working space. 
//     The final elements of BMAT and ZMAT are set in a well-conditioned way
//       to the values that are appropriate for the new interpolation points.
//     The elements of GOPT, HQ and PQ are also revised to the values that are
//       appropriate to the final quadratic model.
//
//     Set some constants.
//
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      NP=N+1
      SFRAC=HALF/DFLOAT(NP)
      NPTM=NPT-NP
//
//     Shift the interpolation points so that XOPT becomes the origin, and set
//     the elements of ZMAT to zero. The value of SUMPQ is required in the
//     updating of HQ below. The squares of the distances from XOPT to the
//     other interpolation points are set at the end of W. Increments of WINC
//     may be added later to these squares to balance the consideration of
//     the choice of point that is going to become current.
//
      SUMPQ=ZERO
      WINC=ZERO
      DO 20 K=1,NPT
      DISTSQ=ZERO
      DO 10 J=1,N
      XPT(K,J)=XPT(K,J)-XOPT[J]
   10 DISTSQ=DISTSQ+XPT(K,J)**2
      SUMPQ=SUMPQ+PQ[K]
      W(NDIM+K)=DISTSQ
      WINC=Math.Max(WINC,DISTSQ)
      DO 20 J=1,NPTM
   20 ZMAT(K,J)=ZERO
//
//     Update HQ so that HQ and PQ define the second derivatives of the model
//     after XBASE has been shifted to the trust region centre.
//
      IH=0
      DO 40 J=1,N
      W[J]=HALF*SUMPQ*XOPT[J]
      DO 30 K=1,NPT
   30 W[J]=W[J]+PQ[K]*XPT(K,J)
      DO 40 I=1,J
      IH=IH+1
   40 HQ(IH)=HQ(IH)+W[I]*XOPT[J]+W[J]*XOPT[I]
//
//     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
//     also set the elements of PTSAUX.
//
      DO 50 J=1,N
      XBASE[J]=XBASE[J]+XOPT[J]
      SL[J]=SL[J]-XOPT[J]
      SU[J]=SU[J]-XOPT[J]
      XOPT[J]=ZERO
      PTSAUX(1,J)=Math.Min(DELTA,SU[J])
      PTSAUX(2,J)=Math.Max(-DELTA,SL[J])
      if (PTSAUX(1,J)+PTSAUX(2,J) < ZERO) {
          TEMP=PTSAUX(1,J)
          PTSAUX(1,J)=PTSAUX(2,J)
          PTSAUX(2,J)=TEMP
      }
      if (DABS(PTSAUX(2,J)) < HALF*DABS(PTSAUX(1,J))) {
          PTSAUX(2,J)=HALF*PTSAUX(1,J)
      }
      DO 50 I=1,NDIM
   50 BMAT(I,J)=ZERO
      FBASE=FVAL(KOPT)
//
//     Set the identifiers of the artificial interpolation points that are
//     along a coordinate direction from XOPT, and set the corresponding
//     nonzero elements of BMAT and ZMAT.
//
      PTSID(1)=SFRAC
      DO 60 J=1,N
      JP=J+1
      JPN=JP+N
      PTSID(JP)=DFLOAT[J]+SFRAC
      if (JPN <= NPT) {
          PTSID(JPN)=DFLOAT[J]/DFLOAT(NP)+SFRAC
          TEMP=ONE/(PTSAUX(1,J)-PTSAUX(2,J))
          BMAT(JP,J)=-TEMP+ONE/PTSAUX(1,J)
          BMAT(JPN,J)=TEMP+ONE/PTSAUX(2,J)
          BMAT(1,J)=-BMAT(JP,J)-BMAT(JPN,J)
          ZMAT(1,J)=Math.Sqrt(2.0D0)/DABS(PTSAUX(1,J)*PTSAUX(2,J))
          ZMAT(JP,J)=ZMAT(1,J)*PTSAUX(2,J)*TEMP
          ZMAT(JPN,J)=-ZMAT(1,J)*PTSAUX(1,J)*TEMP
      ELSE
          BMAT(1,J)=-ONE/PTSAUX(1,J)
          BMAT(JP,J)=ONE/PTSAUX(1,J)
          BMAT(J+NPT,J)=-HALF*PTSAUX(1,J)**2
      }
   60 CONTINUE
//
//     Set any remaining identifiers with their nonzero elements of ZMAT.
//
      if (NPT >= N+NP) {
          DO 70 K=2*NP,NPT
          IW=(DFLOAT(K-NP)-HALF)/DFLOAT(N)
          IP=K-NP-IW*N
          IQ=IP+IW
          if (IQ > N) IQ=IQ-N
          PTSID[K]=DFLOAT(IP)+DFLOAT(IQ)/DFLOAT(NP)+SFRAC
          TEMP=ONE/(PTSAUX(1,IP)*PTSAUX(1,IQ))
          ZMAT(1,K-NP)=TEMP
          ZMAT(IP+1,K-NP)=-TEMP
          ZMAT(IQ+1,K-NP)=-TEMP
   70     ZMAT(K,K-NP)=TEMP
      }
      NREM=NPT
      KOLD=1
      KNEW=KOPT
//
//     Reorder the provisional points in the way that exchanges PTSID(KOLD)
//     with PTSID(KNEW).
//
   80 DO 90 J=1,N
      TEMP=BMAT(KOLD,J)
      BMAT(KOLD,J)=BMAT(KNEW,J)
   90 BMAT(KNEW,J)=TEMP
      DO 100 J=1,NPTM
      TEMP=ZMAT(KOLD,J)
      ZMAT(KOLD,J)=ZMAT(KNEW,J)
  100 ZMAT(KNEW,J)=TEMP
      PTSID(KOLD)=PTSID(KNEW)
      PTSID(KNEW)=ZERO
      W(NDIM+KNEW)=ZERO
      NREM=NREM-1
      if (KNEW != KOPT) {
          TEMP=VLAG(KOLD)
          VLAG(KOLD)=VLAG(KNEW)
          VLAG(KNEW)=TEMP
//
//     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
//     interpolation point can be changed from provisional to original. The
//     branch to label 350 occurs if all the original points are reinstated.
//     The nonnegative values of W(NDIM+K) are required in the search below.
//
          CALL UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
          if (NREM == 0) goto 350
          DO 110 K=1,NPT
  110     W(NDIM+K)=DABS(W(NDIM+K))
      }
//
//     Pick the index KNEW of an original interpolation point that has not
//     yet replaced one of the provisional interpolation points, giving
//     attention to the closeness to XOPT and to previous tries with KNEW.
//
  120 DSQMIN=ZERO
      DO 130 K=1,NPT
      if (W(NDIM+K) > ZERO) {
          if (DSQMIN == ZERO || W(NDIM+K) < DSQMIN) {
              KNEW=K
              DSQMIN=W(NDIM+K)
          }
      }
  130 CONTINUE
      if (DSQMIN == ZERO) goto 260
//
//     Form the W-vector of the chosen original interpolation point.
//
      DO 140 J=1,N
  140 W(NPT+J)=XPT(KNEW,J)
      DO 160 K=1,NPT
      SUM=ZERO
      if (K == KOPT) {
          CONTINUE
      ELSE if (PTSID[K] == ZERO) {
          DO 150 J=1,N
  150     SUM=SUM+W(NPT+J)*XPT(K,J)
      ELSE
          IP=PTSID[K]
          if (IP > 0) SUM=W(NPT+IP)*PTSAUX(1,IP)
          IQ=DFLOAT(NP)*PTSID[K]-DFLOAT(IP*NP)
          if (IQ > 0) {
              IW=1
              if (IP == 0) IW=2
              SUM=SUM+W(NPT+IQ)*PTSAUX(IW,IQ)
          }
      }
  160 W[K]=HALF*SUM*SUM
//
//     Calculate VLAG and BETA for the required updating of the H matrix if
//     XPT(KNEW,.) is reinstated in the set of interpolation points.
//
      DO 180 K=1,NPT
      SUM=ZERO
      DO 170 J=1,N
  170 SUM=SUM+BMAT(K,J)*W(NPT+J)
  180 VLAG[K]=SUM
      BETA=ZERO
      DO 200 J=1,NPTM
      SUM=ZERO
      DO 190 K=1,NPT
  190 SUM=SUM+ZMAT(K,J)*W[K]
      BETA=BETA-SUM*SUM
      DO 200 K=1,NPT
  200 VLAG[K]=VLAG[K]+SUM*ZMAT(K,J)
      BSUM=ZERO
      DISTSQ=ZERO
      DO 230 J=1,N
      SUM=ZERO
      DO 210 K=1,NPT
  210 SUM=SUM+BMAT(K,J)*W[K]
      JP=J+NPT
      BSUM=BSUM+SUM*W(JP)
      DO 220 IP=NPT+1,NDIM
  220 SUM=SUM+BMAT(IP,J)*W(IP)
      BSUM=BSUM+SUM*W(JP)
      VLAG(JP)=SUM
  230 DISTSQ=DISTSQ+XPT(KNEW,J)**2
      BETA=HALF*DISTSQ*DISTSQ+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
//
//     KOLD is set to the index of the provisional interpolation point that is
//     going to be deleted to make way for the KNEW-th original interpolation
//     point. The choice of KOLD is governed by the avoidance of a small value
//     of the denominator in the updating calculation of UPDATE.
//
      DENOM=ZERO
      VLMXSQ=ZERO
      DO 250 K=1,NPT
      if (PTSID[K] != ZERO) {
          HDIAG=ZERO
          DO 240 J=1,NPTM
  240     HDIAG=HDIAG+ZMAT(K,J)**2
          DEN=BETA*HDIAG+VLAG[K]**2
          if (DEN > DENOM) {
              KOLD=K
              DENOM=DEN
          }
      }
  250 VLMXSQ=Math.Max(VLMXSQ,VLAG[K]**2)
      if (DENOM <= 1.0D-2*VLMXSQ) {
          W(NDIM+KNEW)=-W(NDIM+KNEW)-WINC
          goto 120
      }
      goto 80
//
//     When label 260 is reached, all the final positions of the interpolation
//     points have been chosen although any changes have not been included yet
//     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
//     from the shift of XBASE, the updating of the quadratic model remains to
//     be done. The following cycle through the new interpolation points begins
//     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
//     except that a RETURN occurs if MAXFUN prohibits another value of F.
//
  260 DO 340 KPT=1,NPT
      if (PTSID(KPT) == ZERO) goto 340
      if (NF >= MAXFUN) {
          NF=-1
          goto 350
      }
      IH=0
      DO 270 J=1,N
      W[J]=XPT(KPT,J)
      XPT(KPT,J)=ZERO
      TEMP=PQ(KPT)*W[J]
      DO 270 I=1,J
      IH=IH+1
  270 HQ(IH)=HQ(IH)+TEMP*W[I]
      PQ(KPT)=ZERO
      IP=PTSID(KPT)
      IQ=DFLOAT(NP)*PTSID(KPT)-DFLOAT(IP*NP)
      if (IP > 0) {
          XP=PTSAUX(1,IP)
          XPT(KPT,IP)=XP
      }
      if (IQ > 0) {
          XQ=PTSAUX(1,IQ)
          if (IP == 0) XQ=PTSAUX(2,IQ)
          XPT(KPT,IQ)=XQ
      }
//
//     Set VQUAD to the value of the current model at the new point.
//
      VQUAD=FBASE
      if (IP > 0) {
          IHP=(IP+IP*IP)/2
          VQUAD=VQUAD+XP*(GOPT(IP)+HALF*XP*HQ(IHP))
      }
      if (IQ > 0) {
          IHQ=(IQ+IQ*IQ)/2
          VQUAD=VQUAD+XQ*(GOPT(IQ)+HALF*XQ*HQ(IHQ))
          if (IP > 0) {
              IW=MAX0(IHP,IHQ)-IABS(IP-IQ)
              VQUAD=VQUAD+XP*XQ*HQ(IW)
          }
      }
      DO 280 K=1,NPT
      TEMP=ZERO
      if (IP > 0) TEMP=TEMP+XP*XPT(K,IP)
      if (IQ > 0) TEMP=TEMP+XQ*XPT(K,IQ)
  280 VQUAD=VQUAD+HALF*PQ[K]*TEMP*TEMP
//
//     Calculate F at the new interpolation point, and set DIFF to the factor
//     that is going to multiply the KPT-th Lagrange function when the model
//     is updated to provide interpolation to the new function value.
//
      DO 290 I=1,N
      W[I]=Math.Min(Math.Max(XL[I],XBASE[I]+XPT(KPT,I)),XU[I])
      if (XPT(KPT,I) == SL[I]) W[I]=XL[I]
      if (XPT(KPT,I) == SU[I]) W[I]=XU[I]
  290 CONTINUE
      NF=NF+1
      CALL CALFUN (N,W,F)
      if (IPRINT == 3) {
          PRINT 300, NF,F,(W[I],I=1,N)
  300     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1      '    The corresponding X is:'/(2X,5D15.6))
      }
      FVAL(KPT)=F
      if (F < FVAL(KOPT)) KOPT=KPT
      DIFF=F-VQUAD
//
//     Update the quadratic model. The RETURN from the subroutine occurs when
//     all the new interpolation points are included in the model.
//
      DO 310 I=1,N
  310 GOPT[I]=GOPT[I]+DIFF*BMAT(KPT,I)
      DO 330 K=1,NPT
      SUM=ZERO
      DO 320 J=1,NPTM
  320 SUM=SUM+ZMAT(K,J)*ZMAT(KPT,J)
      TEMP=DIFF*SUM
      if (PTSID[K] == ZERO) {
          PQ[K]=PQ[K]+TEMP
      ELSE
          IP=PTSID[K]
          IQ=DFLOAT(NP)*PTSID[K]-DFLOAT(IP*NP)
          IHQ=(IQ*IQ+IQ)/2
          if (IP == 0) {
              HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(2,IQ)**2
          ELSE
              IHP=(IP*IP+IP)/2
              HQ(IHP)=HQ(IHP)+TEMP*PTSAUX(1,IP)**2
              if (IQ > 0) {
                  HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(1,IQ)**2
                  IW=MAX0(IHP,IHQ)-IABS(IQ-IP)
                  HQ(IW)=HQ(IW)+TEMP*PTSAUX(1,IP)*PTSAUX(1,IQ)
              }
          }
      }
  330 CONTINUE
      PTSID(KPT)=ZERO
  340 CONTINUE
  350 RETURN
      END


      SUBROUTINE TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,
     1  XNEW,D,GNEW,XBDI,S,HS,HRED,DSQ,CRVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XPT(NPT,*),XOPT(*),GOPT(*),HQ(*),PQ(*),SL(*),SU(*),
     1  XNEW(*),D(*),GNEW(*),XBDI(*),S(*),HS(*),HRED(*)
//
//     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
//       meanings as the corresponding arguments of BOBYQB.
//     DELTA is the trust region radius for the present calculation, which
//       seeks a small value of the quadratic model within distance DELTA of
//       XOPT subject to the bounds on the variables.
//     XNEW will be set to a new vector of variables that is approximately
//       the one that minimizes the quadratic model within the trust region
//       subject to the SL and SU constraints on the variables. It satisfies
//       as equations the bounds that become active during the calculation.
//     D is the calculated trial step from XOPT, generated iteratively from an
//       initial value of zero. Thus XNEW is XOPT+D after the final iteration.
//     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
//       when D is updated.
//     XBDI is a working space vector. For I=1,2,...,N, the element XBDI[I] is
//       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
//       I-th variable has become fixed at a bound, the bound being SL[I] or
//       SU[I] in the case XBDI[I]=-1.0 or XBDI[I]=1.0, respectively. This
//       information is accumulated during the construction of XNEW.
//     The arrays S, HS and HRED are also used for working space. They hold the
//       current search direction, and the changes in the gradient of Q along S
//       and the reduced D, respectively, where the reduced D is the same as D,
//       except that the components of the fixed variables are zero.
//     DSQ will be set to the square of the length of XNEW-XOPT.
//     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise
//       it is set to the least curvature of H that occurs in the conjugate
//       gradient searches that are not restricted by any constraints. The
//       value CRVMIN=-1.0D0 is set, however, if all of these searches are
//       constrained.
//
//     A version of the truncated conjugate gradient is applied. If a line
//     search is restricted by a constraint, then the procedure is restarted,
//     the values of the variables that are at their bounds being fixed. If
//     the trust region boundary is reached, then further changes may be made
//     to D, each one being in the two dimensional space that is spanned
//     by the current D and the gradient of Q at XOPT+D, staying on the trust
//     region boundary. Termination occurs when the reduction in Q seems to
//     be close to the greatest reduction that can be achieved.
//
//     Set some constants.
//
      HALF=0.5D0
      ONE=1.0D0
      ONEMIN=-1.0D0
      ZERO=0.0D0
//
//     The sign of GOPT[I] gives the sign of the change to the I-th variable
//     that will reduce Q from its value at XOPT. Thus XBDI[I] shows whether
//     or not to fix the I-th variable at one of its bounds initially, with
//     NACT being set to the number of fixed variables. D and GNEW are also
//     set for the first iteration. DELSQ is the upper bound on the sum of
//     squares of the free variables. QRED is the reduction in Q so far.
//
      ITERC=0
      NACT=0
      SQSTP=ZERO
      DO 10 I=1,N
      XBDI[I]=ZERO
      if (XOPT[I] <= SL[I]) {
          if (GOPT[I] >= ZERO) XBDI[I]=ONEMIN
      ELSE if (XOPT[I] >= SU[I]) {
          if (GOPT[I] <= ZERO) XBDI[I]=ONE
      }
      if (XBDI[I] != ZERO) NACT=NACT+1
      D[I]=ZERO
   10 GNEW[I]=GOPT[I]
      DELSQ=DELTA*DELTA
      QRED=ZERO
      CRVMIN=ONEMIN
//
//     Set the next search direction of the conjugate gradient method. It is
//     the steepest descent direction initially and when the iterations are
//     restarted because a variable has just been fixed by a bound, and of
//     course the components of the fixed variables are zero. ITERMAX is an
//     upper bound on the indices of the conjugate gradient iterations.
//
   20 BETA=ZERO
   30 STEPSQ=ZERO
      DO 40 I=1,N
      if (XBDI[I] != ZERO) {
          S[I]=ZERO
      ELSE if (BETA == ZERO) {
          S[I]=-GNEW[I]
      ELSE
          S[I]=BETA*S[I]-GNEW[I]
      }
   40 STEPSQ=STEPSQ+S[I]**2
      if (STEPSQ == ZERO) goto 190
      if (BETA == ZERO) {
          GREDSQ=STEPSQ
          ITERMAX=ITERC+N-NACT
      }
      if (GREDSQ*DELSQ <= 1.0D-4*QRED*QRED) GO TO 190
//
//     Multiply the search direction by the second derivative matrix of Q and
//     calculate some scalars for the choice of steplength. Then set BLEN to
//     the length of the the step to the trust region boundary and STPLEN to
//     the steplength, ignoring the simple bounds.
//
      goto 210
   50 RESID=DELSQ
      DS=ZERO
      SHS=ZERO
      DO 60 I=1,N
      if (XBDI[I] == ZERO) {
          RESID=RESID-D[I]**2
          DS=DS+S[I]*D[I]
          SHS=SHS+S[I]*HS[I]
      }
   60 CONTINUE
      if (RESID <= ZERO) goto 90
      TEMP=Math.Sqrt(STEPSQ*RESID+DS*DS)
      if (DS < ZERO) {
          BLEN=(TEMP-DS)/STEPSQ
      ELSE
          BLEN=RESID/(TEMP+DS)
      }
      STPLEN=BLEN
      if (SHS > ZERO) {
          STPLEN=Math.Min(BLEN,GREDSQ/SHS)
      }
      
//
//     Reduce STPLEN if necessary in order to preserve the simple bounds,
//     letting IACT be the index of the new constrained variable.
//
      IACT=0
      DO 70 I=1,N
      if (S[I] != ZERO) {
          XSUM=XOPT[I]+D[I]
          if (S[I] > ZERO) {
              TEMP=(SU[I]-XSUM)/S[I]
          ELSE
              TEMP=(SL[I]-XSUM)/S[I]
          }
          if (TEMP < STPLEN) {
              STPLEN=TEMP
              IACT=I
          }
      }
   70 CONTINUE
//
//     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
//
      SDEC=ZERO
      if (STPLEN > ZERO) {
          ITERC=ITERC+1
          TEMP=SHS/STEPSQ
          if (IACT == 0 && TEMP > ZERO) {
              CRVMIN=Math.Min(CRVMIN,TEMP)
              if (CRVMIN == ONEMIN) CRVMIN=TEMP
          } 
          GGSAV=GREDSQ
          GREDSQ=ZERO
          DO 80 I=1,N
          GNEW[I]=GNEW[I]+STPLEN*HS[I]
          if (XBDI[I] == ZERO) GREDSQ=GREDSQ+GNEW[I]**2
   80     D[I]=D[I]+STPLEN*S[I]
          SDEC=Math.Max(STPLEN*(GGSAV-HALF*STPLEN*SHS),ZERO)
          QRED=QRED+SDEC
      }
//
//     Restart the conjugate gradient method if it has hit a new bound.
//
      if (IACT > 0) {
          NACT=NACT+1
          XBDI(IACT)=ONE
          if (S(IACT) < ZERO) XBDI(IACT)=ONEMIN
          DELSQ=DELSQ-D(IACT)**2
          if (DELSQ <= ZERO) goto 90
          goto 20
      }
//
//     If STPLEN is less than BLEN, then either apply another conjugate
//     gradient iteration or RETURN.
//
      if (STPLEN < BLEN) {
          if (ITERC == ITERMAX) goto 190
          if (SDEC <= 0.01D0*QRED) goto 190
          BETA=GREDSQ/GGSAV
          goto 30
      }
   90 CRVMIN=ZERO
//
//     Prepare for the alternative iteration by calculating some scalars
//     and by multiplying the reduced D by the second derivative matrix of
//     Q, where S holds the reduced D in the call of GGMULT.
//
  100 if (NACT >= N-1) goto 190
      DREDSQ=ZERO
      DREDG=ZERO
      GREDSQ=ZERO
      DO 110 I=1,N
      if (XBDI[I] == ZERO) {
          DREDSQ=DREDSQ+D[I]**2
          DREDG=DREDG+D[I]*GNEW[I]
          GREDSQ=GREDSQ+GNEW[I]**2
          S[I]=D[I]
      ELSE
          S[I]=ZERO
      }
  110 CONTINUE
      ITCSAV=ITERC
      goto 210
//
//     Let the search direction S be a linear combination of the reduced D
//     and the reduced G that is orthogonal to the reduced D.
//
  120 ITERC=ITERC+1
      TEMP=GREDSQ*DREDSQ-DREDG*DREDG
      if (TEMP <= 1.0D-4*QRED*QRED) goto 190
      TEMP=Math.Sqrt(TEMP)
      DO 130 I=1,N
      if (XBDI[I] == ZERO) {
          S[I]=(DREDG*D[I]-DREDSQ*GNEW[I])/TEMP
      ELSE
          S[I]=ZERO
      }
  130 CONTINUE
      SREDG=-TEMP
//
//     By considering the simple bounds on the variables, calculate an upper
//     bound on the tangent of half the angle of the alternative iteration,
//     namely ANGBD, except that, if already a free variable has reached a
//     bound, there is a branch back to label 100 after fixing that variable.
//
      ANGBD=ONE
      IACT=0
      DO 140 I=1,N
      if (XBDI[I] == ZERO) {
          TEMPA=XOPT[I]+D[I]-SL[I]
          TEMPB=SU[I]-XOPT[I]-D[I]
          if (TEMPA <= ZERO) {
              NACT=NACT+1
              XBDI[I]=ONEMIN
              goto 100
          ELSE if (TEMPB <= ZERO) {
              NACT=NACT+1
              XBDI[I]=ONE
              goto 100
          }
          RATIO=ONE
          SSQ=D[I]**2+S[I]**2
          TEMP=SSQ-(XOPT[I]-SL[I])**2
          if (TEMP > ZERO) {
              TEMP=Math.Sqrt(TEMP)-S[I]
              if (ANGBD*TEMP > TEMPA) {
                  ANGBD=TEMPA/TEMP
                  IACT=I
                  XSAV=ONEMIN
              }
          }
          TEMP=SSQ-(SU[I]-XOPT[I])**2
          if (TEMP > ZERO) {
              TEMP=Math.Sqrt(TEMP)+S[I]
              if (ANGBD*TEMP > TEMPB) {
                  ANGBD=TEMPB/TEMP
                  IACT=I
                  XSAV=ONE
              }
          }
      }
  140 CONTINUE
//
//     Calculate HHD and some curvatures for the alternative iteration.
//
      goto 210
  150 SHS=ZERO
      DHS=ZERO
      DHD=ZERO
      DO 160 I=1,N
      if (XBDI[I] == ZERO) {
          SHS=SHS+S[I]*HS[I]
          DHS=DHS+D[I]*HS[I]
          DHD=DHD+D[I]*HRED[I]
      }
  160 CONTINUE
//
//     Seek the greatest reduction in Q for a range of equally spaced values
//     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
//     the alternative iteration.
//
      REDMAX=ZERO
      ISAV=0
      REDSAV=ZERO
      IU=17.0D0*ANGBD+3.1D0
      DO 170 I=1,IU
      ANGT=ANGBD*DFLOAT[I]/DFLOAT(IU)
      STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
      TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
      REDNEW=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
      if (REDNEW > REDMAX) {
          REDMAX=REDNEW
          ISAV=I
          RDPREV=REDSAV
      ELSE if (I == ISAV+1) {
          RDNEXT=REDNEW
      }
  170 REDSAV=REDNEW
//
//     Return if the reduction is zero. Otherwise, set the sine and cosine
//     of the angle of the alternative iteration, and calculate SDEC.
//
      if (ISAV == 0) goto 190
      if (ISAV < IU) {
          TEMP=(RDNEXT-RDPREV)/(REDMAX+REDMAX-RDPREV-RDNEXT)
          ANGT=ANGBD*(DFLOAT(ISAV)+HALF*TEMP)/DFLOAT(IU)
      }
      CTH=(ONE-ANGT*ANGT)/(ONE+ANGT*ANGT)
      STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
      TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
      SDEC=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
      if (SDEC <= ZERO) goto 190
//
//     Update GNEW, D and HRED. If the angle of the alternative iteration
//     is restricted by a bound on a free variable, that variable is fixed
//     at the bound.
//
      DREDG=ZERO
      GREDSQ=ZERO
      DO 180 I=1,N
      GNEW[I]=GNEW[I]+(CTH-ONE)*HRED[I]+STH*HS[I]
      if (XBDI[I] == ZERO) {
          D[I]=CTH*D[I]+STH*S[I]
          DREDG=DREDG+D[I]*GNEW[I]
          GREDSQ=GREDSQ+GNEW[I]**2
      }
  180 HRED[I]=CTH*HRED[I]+STH*HS[I]
      QRED=QRED+SDEC
      if (IACT > 0 && ISAV == IU) {
          NACT=NACT+1
          XBDI(IACT)=XSAV
          goto 100
      }
//
//     If SDEC is sufficiently small, then RETURN after setting XNEW to
//     XOPT+D, giving careful attention to the bounds.
//
      if (SDEC > 0.01D0*QRED) goto 120
  190 DSQ=ZERO
      DO 200 I=1,N
      XNEW[I]=Math.Max(Math.Min(XOPT[I]+D[I],SU[I]),SL[I])
      if (XBDI[I] == ONEMIN) XNEW[I]=SL[I]
      if (XBDI[I] == ONE) XNEW[I]=SU[I]
      D[I]=XNEW[I]-XOPT[I]
  200 DSQ=DSQ+D[I]**2
      RETURN
 
//     The following instructions multiply the current S-vector by the second
//     derivative matrix of the quadratic model, putting the product in HS.
//     They are reached from three different parts of the software above and
//     they can be regarded as an external subroutine.
//
  210 IH=0
      DO 220 J=1,N
      HS[J]=ZERO
      DO 220 I=1,J
      IH=IH+1
      if (I < J) HS[J]=HS[J]+HQ(IH)*S[I]
  220 HS[I]=HS[I]+HQ(IH)*S[J]
      DO 250 K=1,NPT
      if (PQ[K] != ZERO) {
          TEMP=ZERO
          DO 230 J=1,N
  230     TEMP=TEMP+XPT(K,J)*S[J]
          TEMP=TEMP*PQ[K]
          DO 240 I=1,N
  240     HS[I]=HS[I]+TEMP*XPT(K,I)
      }
  250 CONTINUE
      if (CRVMIN != ZERO) goto 50
      if (ITERC > ITCSAV) goto 150
      DO 260 I=1,N
  260 HRED[I]=HS[I]
      goto 120
      END


      SUBROUTINE UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,
     1  KNEW,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)
//
//     The arrays BMAT and ZMAT are updated, as required by the new position
//     of the interpolation point that has the index KNEW. The vector VLAG has
//     N+NPT components, set on entry to the first NPT and last N components
//     of the product Hw in equation (4.11) of the Powell (2006) paper on
//     NEWUOA. Further, BETA is set on entry to the value of the parameter
//     with that name, and DENOM is set to the denominator of the updating
//     formula. Elements of ZMAT may be treated as zero if their moduli are
//     at most ZTEST. The first NDIM elements of W are used for working space.
//
//     Set some constants.
//
      ONE=1.0D0
      ZERO=0.0D0
      NPTM=NPT-N-1
      ZTEST=ZERO
      DO 10 K=1,NPT
      DO 10 J=1,NPTM
   10 ZTEST=Math.Max(ZTEST,DABS(ZMAT(K,J)))
      ZTEST=1.0D-20*ZTEST
//
//     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
//
      JL=1
      DO 30 J=2,NPTM
      if (DABS(ZMAT(KNEW,J)) > ZTEST) {
          TEMP=Math.Sqrt(ZMAT(KNEW,1)**2+ZMAT(KNEW,J)**2)
          TEMPA=ZMAT(KNEW,1)/TEMP
          TEMPB=ZMAT(KNEW,J)/TEMP
          DO 20 I=1,NPT
          TEMP=TEMPA*ZMAT(I,1)+TEMPB*ZMAT(I,J)
          ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,1)
   20     ZMAT(I,1)=TEMP
      }
      ZMAT(KNEW,J)=ZERO
   30 CONTINUE
//
//     Put the first NPT components of the KNEW-th column of HLAG into W,
//     and calculate the parameters of the updating formula.
//
      DO 40 I=1,NPT
      W[I]=ZMAT(KNEW,1)*ZMAT(I,1)
   40 CONTINUE
      ALPHA=W(KNEW)
      TAU=VLAG(KNEW)
      VLAG(KNEW)=VLAG(KNEW)-ONE
//
//     Complete the updating of ZMAT.
//
      TEMP=Math.Sqrt(DENOM)
      TEMPB=ZMAT(KNEW,1)/TEMP
      TEMPA=TAU/TEMP
      DO 50 I=1,NPT
   50 ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG[I]
//
//     Finally, update the matrix BMAT.
//
      DO 60 J=1,N
      JP=NPT+J
      W(JP)=BMAT(KNEW,J)
      TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
      TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
      DO 60 I=1,JP
      BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG[I]+TEMPB*W[I]
      if (I > NPT) BMAT(JP,I-NPT)=BMAT(I,J)
   60 CONTINUE
      RETURN
      END
*/