/*
 * velschemers.cpp
 *
 *     Project: Finite Volume Navier-Stokes Solver
 *      Author: Syed Ahmad Raza
 */

// #include <cmath>        // math functions
// #include <iostream>     // functions for input and output to console

#include "constants.h"
#include "velschemers.h"


using namespace Numeric_lib;

void quick(int i, int j, int k)
{
    if (U(i,j,k) > 0)
    {
        ue  = 0.5*(U(i,j,k)     + U(i+1,j,k)  ) -0.125*Dx(i+1)*Dx(i+1)/Dxs(i)
              *(  (U(i+1,j,k)   - U(i,j,k)    ) / Dx(i+1)
                - (U(i,j,k)     - U(i-1,j,k)  ) / Dx(i)   );

        uw  = 0.5*(U(i-1,j,k)   + U(i,j,k)    ) -0.125*Dx(i)*Dx(i)/Dxs(i-1)
              *(  (U(i,j,k)     - U(i-1,j,k)  ) / Dx(i)
                - (U(i-1,j,k)   - U(i-2,j,k)  ) / Dx(i-1) );

        ve  = 0.5*(V(i,j,k)     + V(i+1,j,k)  ) -0.125*Dxs(i)*Dxs(i)/Dx(i)
              *(  (V(i+1,j,k)   - V(i,j,k)    ) / Dxs(i)
                - (V(i,j,k)     - V(i-1,j,k)  ) / Dxs(i-1));

        vw  = 0.5*(V(i-1,j,k)   + V(i,j,k)    ) -0.125*Dxs(i-1)*Dxs(i-1)/Dx(i-1)
              *(  (V(i,j,k)     - V(i-1,j,k)  ) / Dxs(i-1)
                - (V(i-1,j,k)   - V(i-2,j,k)  ) / Dxs(i-2));

        we  = 0.5*(W(i,j,k)     + W(i+1,j,k)  ) -0.125*Dxs(i)*Dxs(i)/Dx(i)
              *(  (W(i+1,j,k)   - W(i,j,k)    ) / Dxs(i)
                - (W(i,j,k)     - W(i-1,j,k)  ) / Dxs(i-1));

        ww  = 0.5*(W(i-1,j,k)   + W(i,j,k)    ) -0.125*Dxs(i-1)*Dxs(i-1)/Dx(i-1)
              *(  (W(i,j,k)     - W(i-1,j,k)  ) / Dxs(i-1)
                - (W(i-1,j,k)   - W(i-2,j,k)  ) / Dxs(i-2));

        vnu = 0.5*(V(i,j,k)     + V(i+1,j,k)  ) -0.125*Dxs(i)*Dxs(i)/Dx(i)
              *(  (V(i+1,j,k)   - V(i,j,k)    ) / Dxs(i)
                - (V(i,j,k)     - V(i-1,j,k)  ) / Dxs(i-1));

        vsu = 0.5*(V(i,j-1,k)   + V(i+1,j-1,k)) -0.125*Dxs(i)*Dxs(i)/Dx(i)
              *(  (V(i+1,j-1,k) - V(i,j-1,k)  ) / Dxs(i)
                - (V(i,j-1,k)   - V(i-1,j-1,k)) / Dxs(i-1));

        wfu = 0.5*(W(i,j,k)     + W(i+1,j,k)  ) -0.125*Dxs(i)*Dxs(i)/Dx(i)
              *(  (W(i+1,j,k)   - W(i,j,k)    ) / Dxs(i)
                - (W(i,j,k)     - W(i-1,j,k)  ) / Dxs(i-1));

        wbu = 0.5*(W(i,j,k-1)   + W(i+1,j,k-1)) -0.125*Dxs(i)*Dxs(i)/Dx(i)
              *(  (W(i+1,j,k-1) - W(i,j,k-1)  ) / Dxs(i)
                - (W(i,j,k-1)   - W(i-1,j,k-1)) / Dxs(i-1));
    }
    else
    {
        ue  = 0.5*(U(i,j,k)     + U(i+1,j,k)  ) -0.125*Dx(i+1)*Dx(i+1)/Dxs(i+1)
              *(  (U(i+2,j,k)   - U(i+1,j,k)  ) / Dx(i+2)
                - (U(i+1,j,k)   - U(i,j,k)    ) / Dx(i+1) );

        uw  = 0.5*(U(i-1,j,k)   + U(i,j,k)    ) -0.125*Dx(i)*Dx(i)/Dxs(i)
              *(  (U(i+1,j,k)   - U(i,j,k)    ) / Dx(i+1)
                - (U(i,j,k)     - U(i-1,j,k)  ) / Dx(i)   );

        ve  = 0.5*(V(i,j,k)     + V(i+1,j,k)  ) -0.125*Dxs(i)*Dxs(i)/Dx(i+1)
              *(  (V(i+2,j,k)   - V(i+1,j,k)  ) / Dxs(i+1)
                - (V(i+1,j,k)   - V(i,j,k)    ) / Dxs(i)  );

        vw  = 0.5*(V(i-1,j,k)   + V(i,j,k)    ) -0.125*Dxs(i-1)*Dxs(i-1)/Dx(i)
              *(  (V(i+1,j,k)   - V(i,j,k)    ) / Dxs(i)
                - (V(i,j,k)     - V(i-1,j,k)  ) / Dxs(i-1));

        we  = 0.5*(W(i,j,k)     + W(i+1,j,k)  ) -0.125*Dxs(i)*Dxs(i)/Dx(i+1)
              *(  (W(i+2,j,k)   - W(i+1,j,k)  ) / Dxs(i+1)
                - (W(i+1,j,k)   - W(i,j,k)    ) / Dxs(i)  );

        ww  = 0.5*(W(i-1,j,k)   + W(i,j,k)    ) -0.125*Dxs(i-1)*Dxs(i-1)/Dx(i)
              *(  (W(i+1,j,k)   - W(i,j,k)    ) / Dxs(i)
                - (W(i,j,k)     - W(i-1,j,k)  ) / Dxs(i-1));

        vnu = 0.5*(V(i,j,k)     + V(i+1,j,k)  ) -0.125*Dxs(i)*Dxs(i)/Dx(i+1)
              *(  (V(i+2,j,k)   - V(i+1,j,k)  ) / Dxs(i+1)
                - (V(i+1,j,k)   - V(i,j,k)    ) / Dxs(i)  );

        vsu = 0.5*(V(i,j-1,k)   + V(i+1,j-1,k)) -0.125*Dxs(i)*Dxs(i)/Dx(i+1)
              *(  (V(i+2,j-1,k) - V(i+1,j-1,k)) / Dxs(i+1)
                - (V(i+1,j-1,k) - V(i,j-1,k)  ) / Dxs(i)  );

        wfu = 0.5*(W(i,j,k)     + W(i+1,j,k)  ) -0.125*Dxs(i)*Dxs(i)/Dx(i+1)
              *(  (W(i+2,j,k)   - W(i+1,j,k)  ) / Dxs(i+1)
                - (W(i+1,j,k)   - W(i,j,k)    ) / Dxs(i)  );

        wbu = 0.5*(W(i,j,k-1)   + W(i+1,j,k-1)) -0.125*Dxs(i)*Dxs(i)/Dx(i+1)
              *(  (W(i+2,j,k-1) - W(i+1,j,k-1)) / Dxs(i+1)
                - (W(i+1,j,k-1) - W(i,j,k-1)  ) / Dxs(i)  );
    }
    if (V(i,j,k) > 0)
    {
        vn  = 0.5*(V(i,j,k)     + V(i,j+1,k)  ) -0.125*Dy(j+1)*Dy(j+1)/Dys(j)
              *(  (V(i,j+1,k)   - V(i,j,k)    ) / Dy(j+1)
                - (V(i,j,k)     - V(i,j-1,k)  ) / Dy(j)   );

        vs  = 0.5*(V(i,j-1,k)   + V(i,j,k)    ) -0.125*Dy(j)*Dy(j)/Dys(j-1)
              *(  (V(i,j,k)     - V(i,j-1,k)  ) / Dy(j)
                - (V(i,j-1,k)   - V(i,j-2,k)  ) / Dy(j-1) );

        un  = 0.5*(U(i,j,k)     + U(i,j+1,k)  ) -0.125*Dys(j)*Dys(j)/Dy(j)
              *(  (U(i,j+1,k)   - U(i,j,k)    ) / Dys(j)
                - (U(i,j,k)     - U(i,j-1,k)  ) / Dys(j-1));

        us  = 0.5*(U(i,j-1,k)   + U(i,j,k)    ) -0.125*Dys(j-1)*Dys(j-1)/Dy(j-1)
              *(  (U(i,j,k)     - U(i,j-1,k)  ) / Dys(j-1)
                - (U(i,j-1,k)   - U(i,j-2,k)  ) / Dys(j-2));

        wn  = 0.5*(W(i,j,k)     + W(i,j+1,k)  ) -0.125*Dys(j)*Dys(j)/Dy(j)
              *(  (W(i,j+1,k)   - W(i,j,k)    ) / Dys(j)
                - (W(i,j,k)     - W(i,j-1,k)  ) / Dys(j-1));

        ws  = 0.5*(W(i,j-1,k)   + W(i,j,k)    ) -0.125*Dys(j-1)*Dys(j-1)/Dy(j-1)
              *(  (W(i,j,k)     - W(i,j-1,k)  ) / Dys(j-1)
                - (W(i,j-1,k)   - W(i,j-2,k)  ) / Dys(j-2));

        uev = 0.5*(U(i,j,k)     + U(i,j+1,k)  ) -0.125*Dys(j)*Dys(j)/Dy(j)
              *(  (U(i,j+1,k)   - U(i,j,k)    ) / Dys(j)
                - (U(i,j,k)     - U(i,j-1,k)  ) / Dys(j-1));

        uwv = 0.5*(U(i-1,j,k)   + U(i-1,j+1,k)) -0.125*Dys(j)*Dys(j)/Dy(j)
              *(  (U(i-1,j+1,k) - U(i-1,j,k)  ) / Dys(j)
                - (U(i-1,j,k)   - U(i-1,j-1,k)) / Dys(j-1));

        wfv = 0.5*(W(i,j,k)     + W(i,j+1,k)  ) -0.125*Dys(j)*Dys(j)/Dy(j)
              *(  (W(i,j+1,k)   - W(i,j,k)    ) / Dys(j)
                - (W(i,j,k)     - W(i,j-1,k)  ) / Dys(j-1));

        wbv = 0.5*(W(i,j,k-1)   + W(i,j+1,k-1)) -0.125*Dys(j)*Dys(j)/Dy(j)
              *(  (W(i,j+1,k-1) - W(i,j,k-1)  ) / Dys(j)
                - (W(i,j,k-1)   - W(i,j-1,k-1)) / Dys(j-1));
    }
    else
    {
        vn  = 0.5*(V(i,j,k)     + V(i,j+1,k)  ) -0.125*Dy(j+1)*Dy(j+1)/Dys(j+1)
              *(  (V(i,j+2,k)   - V(i,j+1,k)  ) / Dy(j+2)
                - (V(i,j+1,k)   - V(i,j,k)    ) / Dy(j+1) );

        vs  = 0.5*(V(i,j-1,k)   + V(i,j,k)    ) -0.125*Dy(j)*Dy(j)/Dys(j)
              *(  (V(i,j+1,k)   - V(i,j,k)    ) / Dy(j+1)
                - (V(i,j,k)     - V(i,j-1,k)  ) / Dy(j)   );

        un  = 0.5*(U(i,j,k)     + U(i,j+1,k)  ) -0.125*Dys(j)*Dys(j)/Dy(j+1)
              *(  (U(i,j+2,k)   - U(i,j+1,k)  ) / Dys(j+1)
                - (U(i,j+1,k)   - U(i,j,k)    ) / Dys(j)  );

        us  = 0.5*(U(i,j-1,k)   + U(i,j,k)    ) -0.125*Dys(j-1)*Dys(j-1)/Dy(j)
              *(  (U(i,j+1,k)   - U(i,j,k)    ) / Dys(j)
                - (U(i,j,k)     - U(i,j-1,k)  ) / Dys(j-1));

        wn  = 0.5*(W(i,j,k)     + W(i,j+1,k)  ) -0.125*Dys(j)*Dys(j)/Dy(j+1)
              *(  (W(i,j+2,k)   - W(i,j+1,k)  ) / Dys(j+1)
                - (W(i,j+1,k)   - W(i,j,k)    ) / Dys(j)  );

        ws  = 0.5*(W(i,j-1,k)   + W(i,j,k)    ) -0.125*Dys(j-1)*Dys(j-1)/Dy(j)
              *(  (W(i,j+1,k)   - W(i,j,k)    ) / Dys(j)
                - (W(i,j,k)     - W(i,j-1,k)  ) / Dys(j-1));

        uev = 0.5*(U(i,j,k)     + U(i,j+1,k)  ) -0.125*Dys(j)*Dys(j)/Dy(j+1)
              *(  (U(i,j+2,k)   - U(i,j+1,k)  ) / Dys(j+1)
                - (U(i,j+1,k)   - U(i,j,k)    ) / Dys(j)  );

        uwv = 0.5*(U(i-1,j,k)   + U(i-1,j+1,k)) -0.125*Dys(j)*Dys(j)/Dy(j+1)
              *(  (U(i-1,j+2,k) - U(i-1,j+1,k)) / Dys(j+1)
                - (U(i-1,j+1,k) - U(i-1,j,k)  ) / Dys(j)  );

        wfv = 0.5*(W(i,j,k)     + W(i,j+1,k)  ) -0.125*Dys(j)*Dys(j)/Dy(j+1)
              *(  (W(i,j+2,k)   - W(i,j+1,k)  ) / Dys(j+1)
                - (W(i,j+1,k)   - W(i,j,k)    ) / Dys(j)  );

        wbv = 0.5*(W(i,j,k-1)   + W(i,j+1,k-1)) -0.125*Dys(j)*Dys(j)/Dy(j+1)
              *(  (W(i,j+2,k-1) - W(i,j+1,k-1)) / Dys(j+1)
                - (W(i,j+1,k-1) - W(i,j,k-1)  ) / Dys(j)  );
    }
    if (W(i,j,k) > 0)
    {
        wf  = 0.5*(W(i,j,k)     + W(i,j,k+1)  ) -0.125*Dz(k+1)*Dz(k+1)/Dzs(k)
              *(  (W(i,j,k+1)   - W(i,j,k)    ) / Dz(k+1)
                - (W(i,j,k)     - W(i,j,k-1)  ) / Dz(k)   );

        wb  = 0.5*(W(i,j,k-1)   + W(i,j,k)    ) -0.125*Dz(k)*Dz(k)/Dzs(k-1)
              *(  (W(i,j,k)     - W(i,j,k-1)  ) / Dz(k)
                - (W(i,j,k-1)   - W(i,j,k-2)  ) / Dz(k-1) );

        uf  = 0.5*(U(i,j,k)     + U(i,j,k+1)  ) -0.125*Dzs(k)*Dzs(k)/Dz(k)
              *(  (U(i,j,k+1)   - U(i,j,k)    ) / Dzs(k)
                - (U(i,j,k)     - U(i,j,k-1)  ) / Dzs(k-1));

        ub  = 0.5*(U(i,j,k-1)   + U(i,j,k)    ) -0.125*Dzs(k-1)*Dzs(k-1)/Dz(k-1)
              *(  (U(i,j,k)     - U(i,j,k-1)  ) / Dzs(k-1)
                - (U(i,j,k-1)   - U(i,j,k-2)  ) / Dzs(k-2));

        vf  = 0.5*(V(i,j,k)     + V(i,j,k+1)  ) -0.125*Dzs(k)*Dzs(k)/Dz(k)
              *(  (V(i,j,k+1)   - V(i,j,k)    ) / Dzs(k)
                - (V(i,j,k)     - V(i,j,k-1)  ) / Dzs(k-1));

        vb  = 0.5*(V(i,j,k-1)   + V(i,j,k)    ) -0.125*Dzs(k-1)*Dzs(k-1)/Dz(k-1)
              *(  (V(i,j,k)     - V(i,j,k-1)  ) / Dzs(k-1)
                - (V(i,j,k-1)   - V(i,j,k-2)  ) / Dzs(k-2));

        uew = 0.5*(U(i,j,k)     + U(i,j,k+1)  ) -0.125*Dzs(k)*Dzs(k)/Dz(k)
              *(  (U(i,j,k+1)   - U(i,j,k)    ) / Dzs(k)
                - (U(i,j,k)     - U(i,j,k-1)  ) / Dzs(k-1));

        uww = 0.5*(U(i-1,j,k)   + U(i-1,j,k+1)) -0.125*Dzs(k)*Dzs(k)/Dz(k)
              *(  (U(i-1,j,k+1) - U(i-1,j,k)  ) / Dzs(k)
                - (U(i-1,j,k)   - U(i-1,j,k-1)) / Dzs(k-1));

        vnw = 0.5*(V(i,j,k)     + V(i,j,k+1)  ) -0.125*Dzs(k)*Dzs(k)/Dz(k)
              *(  (V(i,j,k+1)   - V(i,j,k)    ) / Dzs(k)
                - (V(i,j,k)     - V(i,j,k-1)  ) / Dzs(k-1));

        vsw = 0.5*(V(i,j-1,k)   + V(i,j-1,k+1)) -0.125*Dzs(k)*Dzs(k)/Dz(k)
              *(  (V(i,j-1,k+1) - V(i,j-1,k)  ) / Dzs(k)
                - (V(i,j-1,k)   - V(i,j-1,k-1)) / Dzs(k-1));
    }
    else
    {
        wf  = 0.5*(W(i,j,k)     + W(i,j,k+1)  ) -0.125*Dz(k+1)*Dz(k+1)/Dzs(k+1)
              *(  (W(i,j,k+2)   - W(i,j,k+1)  ) / Dz(k+2)
                - (W(i,j,k+1)   - W(i,j,k)    ) / Dz(k+1) );

        wb  = 0.5*(W(i,j,k-1)   + W(i,j,k)    ) -0.125*Dz(k)*Dz(k)/Dzs(k)
              *(  (W(i,j,k+1)   - W(i,j,k)    ) / Dz(k+1)
                - (W(i,j,k)     - W(i,j,k-1)  ) / Dz(k)   );

        uf  = 0.5*(U(i,j,k)     + U(i,j,k+1)  ) -0.125*Dzs(k)*Dzs(k)/Dz(k+1)
              *(  (U(i,j,k+2)   - U(i,j,k+1)  ) / Dzs(k+1)
                - (U(i,j,k+1)   - U(i,j,k)    ) / Dzs(k)  );

        ub  = 0.5*(U(i,j,k-1)   + U(i,j,k)    ) -0.125*Dzs(k-1)*Dzs(k-1)/Dz(k)
              *(  (U(i,j,k+1)   - U(i,j,k)    ) / Dzs(k)
                - (U(i,j,k)     - U(i,j,k-1)  ) / Dzs(k-1));

        vf  = 0.5*(V(i,j,k)     + V(i,j,k+1)  ) -0.125*Dzs(k)*Dzs(k)/Dz(k+1)
              *(  (V(i,j,k+2)   - V(i,j,k+1)  ) / Dzs(k+1)
                - (V(i,j,k+1)   - V(i,j,k)    ) / Dzs(k)  );

        vb  = 0.5*(V(i,j,k-1)   + V(i,j,k)    ) -0.125*Dzs(k-1)*Dzs(k-1)/Dz(k)
              *(  (V(i,j,k+1)   - V(i,j,k)    ) / Dzs(k)
                - (V(i,j,k)     - V(i,j,k-1)  ) / Dzs(k-1));

        uew = 0.5*(U(i,j,k)     + U(i,j,k+1)  ) -0.125*Dzs(k)*Dzs(k)/Dz(k+1)
              *(  (U(i,j,k+2)   - U(i,j,k+1)  ) / Dzs(k+1)
                - (U(i,j,k+1)   - U(i,j,k)    ) / Dzs(k)  );

        uww = 0.5*(U(i-1,j,k)   + U(i-1,j,k+1)) -0.125*Dzs(k)*Dzs(k)/Dz(k+1)
              *(  (U(i-1,j,k+2) - U(i-1,j,k+1)) / Dzs(k+1)
                - (U(i-1,j,k+1) - U(i-1,j,k)  ) / Dzs(k)  );

        vnw = 0.5*(V(i,j,k)     + V(i,j,k+1)  ) -0.125*Dzs(k)*Dzs(k)/Dz(k+1)
              *(  (V(i,j,k+2)   - V(i,j,k+1)  ) / Dzs(k+1)
                - (V(i,j,k+1)   - V(i,j,k)    ) / Dzs(k)  );

        vsw = 0.5*(V(i,j-1,k)   + V(i,j-1,k+1)) -0.125*Dzs(k)*Dzs(k)/Dz(k+1)
              *(  (V(i,j-1,k+2) - V(i,j-1,k+1)) / Dzs(k+1)
                - (V(i,j-1,k+1) - V(i,j-1,k)  ) / Dzs(k)  );
    }
}

void upwind(int i, int j, int k)
{
    if (U(i,j,k) > 0)
    {
        ue  = U(i,j,k);
        uw  = U(i-1,j,k);
        vnu = V(i,j,k);
        vsu = V(i,j-1,k);
        wfu = W(i,j,k);
        wbu = W(i,j,k-1);
    }
    else
    {
        ue  = U(i+1,j,k);
        uw  = U(i,j,k);
        vnu = V(i+1,j,k);
        vsu = V(i+1,j-1,k);
        wfu = W(i+1,j,k);
        wbu = W(i+1,j,k-1);
    }
    if (V(i,j,k) > 0)
    {
        uev = U(i,j,k);
        uwv = U(i-1,j,k);
        vn  = V(i,j,k);
        vs  = V(i,j-1,k);
        wfv = W(i,j,k);
        wbv = W(i,j,k-1);
    }
    else
    {
        uev = U(i,j+1,k);
        uwv = U(i-1,j+1,k);
        vn  = V(i,j+1,k);
        vs  = V(i,j,k);
        wfv = W(i,j+1,k);
        wbv = W(i,j+1,k-1);
    }
    if (W(i,j,k) > 0)
    {
        uew = U(i,j,k);
        uww = U(i-1,j,k);
        vnw = V(i,j,k);
        vsw = V(i,j-1,k);
        wf  = W(i,j,k);
        wb  = W(i,j,k-1);
    }
    else
    {
        uew = U(i,j+1,k);
        uww = U(i-1,j,k+1);
        vnw = V(i,j,k+1);
        vsw = V(i,j-1,k+1);
        wf  = W(i,j,k+1);
        wb  = W(i,j,k);
    }
    un = ( U(i,j+1,k) + U(i,j,k)   ) / 2.0;
    us = ( U(i,j,k)   + U(i,j-1,k) ) / 2.0;
    uf = ( U(i,j,k+1) + U(i,j,k)   ) / 2.0;
    ub = ( U(i,j,k)   + U(i,j,k-1) ) / 2.0;

    ve = ( V(i+1,j,k) + V(i,j,k)   ) / 2.0;
    vw = ( V(i,j,k)   + V(i-1,j,k) ) / 2.0;
    vf = ( V(i,j,k+1) + V(i,j,k)   ) / 2.0;
    vb = ( V(i,j,k)   + V(i,j,k-1) ) / 2.0;

    we = ( W(i+1,j,k) + W(i,j,k)   ) / 2.0;
    ww = ( W(i,j,k)   + W(i-1,j,k) ) / 2.0;
    wn = ( W(i,j+1,k) + W(i,j,k)   ) / 2.0;
    ws = ( W(i,j,k)   + W(i,j-1,k) ) / 2.0;
}
