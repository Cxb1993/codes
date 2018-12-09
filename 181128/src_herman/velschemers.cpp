/*
 * velschemers.cpp
 *
 *     Project: Finite Volume Navier-Stokes Solver
 *      Author: Syed Ahmad Raza
 */

#include "constants.h"
#include "velschemers.h"


using namespace Numeric_lib;


void quick(int i, int j,\
    double& ue, double& uw, double& un, double& us,\
    double& vnu, double& vsu,\
    double& ve, double& vw, double& vn, double& vs,\
    double& uev, double& uwv)
{
    if (U(i,j) > 0)
    {
        ue  = 0.5*(U(i,j)     + U(i+1,j)  ) -0.125*Dx(i+1)*Dx(i+1)/Dxs(i)
              *(  (U(i+1,j)   - U(i,j)    ) / Dx(i+1)
                - (U(i,j)     - U(i-1,j)  ) / Dx(i)  );

        uw  = 0.5*(U(i-1,j)   + U(i,j)    ) -0.125*Dx(i)*Dx(i)/Dxs(i-1)
              *(  (U(i,j)     - U(i-1,j)  ) / Dx(i)
                - (U(i-1,j)   - U(i-2,j)  ) / Dx(i-1) );

        ve  = 0.5*(V(i,j)     + V(i+1,j)  ) -0.125*Dxs(i)*Dxs(i)/Dx(i)
              *(  (V(i+1,j)   - V(i,j)    ) / Dxs(i)
                - (V(i,j)     - V(i-1,j)  ) / Dxs(i-1));

        vw  = 0.5*(V(i-1,j)   + V(i,j)    ) -0.125*Dxs(i-1)*Dxs(i-1)/Dx(i-1)
              *(  (V(i,j)     - V(i-1,j)  ) / Dxs(i-1)
                - (V(i-1,j)   - V(i-2,j)  ) / Dxs(i-2));

        vnu = 0.5*(V(i,j)     + V(i+1,j)  ) -0.125*Dxs(i)*Dxs(i)/Dx(i)
              *(  (V(i+1,j)   - V(i,j)    ) / Dxs(i)
                - (V(i,j)     - V(i-1,j)  ) / Dxs(i-1));

        vsu = 0.5*(V(i,j-1)   + V(i+1,j-1)) -0.125*Dxs(i)*Dxs(i)/Dx(i)
              *(  (V(i+1,j-1) - V(i,j-1)  ) / Dxs(i)
                - (V(i,j-1)   - V(i-1,j-1)) / Dxs(i-1));
    }
    else
    {
        ue  = 0.5*(U(i,j)     + U(i+1,j)  ) -0.125*Dx(i+1)*Dx(i+1)/Dxs(i+1)
              *(  (U(i+2,j)   - U(i+1,j)  ) / Dx(i+2)
                - (U(i+1,j)   - U(i,j)    ) / Dx(i+1) );

        uw  = 0.5*(U(i-1,j)   + U(i,j)    ) -0.125*Dx(i)*Dx(i)/Dxs(i)
              *(  (U(i+1,j)   - U(i,j)    ) / Dx(i+1)
                - (U(i,j)     - U(i-1,j)  ) / Dx(i)   );

        ve  = 0.5*(V(i,j)     + V(i+1,j)  ) -0.125*Dxs(i)*Dxs(i)/Dx(i+1)
              *(  (V(i+2,j)   - V(i+1,j)  ) / Dxs(i+1)
                - (V(i+1,j)   - V(i,j)    ) / Dxs(i)  );

        vw  = 0.5*(V(i-1,j)   + V(i,j)    ) -0.125*Dxs(i-1)*Dxs(i-1)/Dx(i)
              *(  (V(i+1,j)   - V(i,j)    ) / Dxs(i)
                - (V(i,j)     - V(i-1,j)  ) / Dxs(i-1));

        vnu = 0.5*(V(i,j)     + V(i+1,j)  ) -0.125*Dxs(i)*Dxs(i)/Dx(i+1)
              *(  (V(i+2,j)   - V(i+1,j)  ) / Dxs(i+1)
                - (V(i+1,j)   - V(i,j)    ) / Dxs(i)  );

        vsu = 0.5*(V(i,j-1)   + V(i+1,j-1)) -0.125*Dxs(i)*Dxs(i)/Dx(i+1)
              *(  (V(i+2,j-1) - V(i+1,j-1)) / Dxs(i+1)
                - (V(i+1,j-1) - V(i,j-1)  ) / Dxs(i)  );
    }
    if (V(i,j) > 0)
    {
        vn  = 0.5*(V(i,j)     + V(i,j+1)  ) -0.125*Dy(j+1)*Dy(j+1)/Dys(j)
              *(  (V(i,j+1)   - V(i,j)    ) / Dy(j+1)
                - (V(i,j)    - V(i,j-1)  ) / Dy(j)   );

        vs  = 0.5*(V(i,j-1)   + V(i,j)    ) -0.125*Dy(j)*Dy(j)/Dys(j-1)
              *(  (V(i,j)     - V(i,j-1)  ) / Dy(j)
                - (V(i,j-1)   - V(i,j-2)  ) / Dy(j-1) );

        un  = 0.5*(U(i,j)     + U(i,j+1)  ) -0.125*Dys(j)*Dys(j)/Dy(j)
              *(  (U(i,j+1)   - U(i,j)    ) / Dys(j)
                - (U(i,j)     - U(i,j-1)  ) / Dys(j-1));

        us  = 0.5*(U(i,j-1)   + U(i,j)    ) -0.125*Dys(j-1)*Dys(j-1)/Dy(j-1)
              *(  (U(i,j)     - U(i,j-1)  ) / Dys(j-1)
                - (U(i,j-1)   - U(i,j-2)  ) / Dys(j-2));

        uev = 0.5*(U(i,j)     + U(i,j+1)  ) -0.125*Dys(j)*Dys(j)/Dy(j)
              *(  (U(i,j+1)   - U(i,j)    ) / Dys(j)
                - (U(i,j)     - U(i,j-1)  ) / Dys(j-1));

        uwv = 0.5*(U(i-1,j)   + U(i-1,j+1)) -0.125*Dys(j)*Dys(j)/Dy(j)
              *(  (U(i-1,j+1) - U(i-1,j)  ) / Dys(j)
                - (U(i-1,j)   - U(i-1,j-1)) / Dys(j-1));
    }
    else
    {
        vn  = 0.5*(V(i,j)     + V(i,j+1)  ) -0.125*Dy(j+1)*Dy(j+1)/Dys(j+1)
              *(  (V(i,j+2)   - V(i,j+1)  ) / Dy(j+2)
                - (V(i,j+1)   - V(i,j)    ) / Dy(j+1) );

        vs  = 0.5*(V(i,j-1)   + V(i,j)    ) -0.125*Dy(j)*Dy(j)/Dys(j)
              *(  (V(i,j+1)   - V(i,j)    ) / Dy(j+1)
                - (V(i,j)     - V(i,j-1)  ) / Dy(j)   );

        un  = 0.5*(U(i,j)     + U(i,j+1)  ) -0.125*Dys(j)*Dys(j)/Dy(j+1)
              *(  (U(i,j+2)   - U(i,j+1)  ) / Dys(j+1)
                - (U(i,j+1)   - U(i,j)    ) / Dys(j)  );

        us  = 0.5*(U(i,j-1)   + U(i,j)    ) -0.125*Dys(j-1)*Dys(j-1)/Dy(j)
              *(  (U(i,j+1)   - U(i,j)    ) / Dys(j)
                - (U(i,j)     - U(i,j-1)  ) / Dys(j-1));

        uev = 0.5*(U(i,j)     + U(i,j+1)  ) -0.125*Dys(j)*Dys(j)/Dy(j+1)
              *(  (U(i,j+2)   - U(i,j+1)  ) / Dys(j+1)
                - (U(i,j+1)   - U(i,j)    ) / Dys(j)  );

        uwv = 0.5*(U(i-1,j)   + U(i-1,j+1)) -0.125*Dys(j)*Dys(j)/Dy(j+1)
              *(  (U(i-1,j+2) - U(i-1,j+1)) / Dys(j+1)
                - (U(i-1,j+1) - U(i-1,j)  ) / Dys(j)  );
    }


}

// void upwind(int i, int j, int k)
// {
//     if (U(i,j,k) > 0)
//     {
//         ue  = U(i,j,k);
//         uw  = U(i-1,j,k);
//         vnu = V(i,j,k);
//         vsu = V(i,j-1,k);
//         wfu = W(i,j,k);
//         wbu = W(i,j,k-1);
//     }
//     else
//     {
//         ue  = U(i+1,j,k);
//         uw  = U(i,j,k);
//         vnu = V(i+1,j,k);
//         vsu = V(i+1,j-1,k);
//         wfu = W(i+1,j,k);
//         wbu = W(i+1,j,k-1);
//     }
//     if (V(i,j,k) > 0)
//     {
//         uev = U(i,j,k);
//         uwv = U(i-1,j,k);
//         vn  = V(i,j,k);
//         vs  = V(i,j-1,k);
//         wfv = W(i,j,k);
//         wbv = W(i,j,k-1);
//     }
//     else
//     {
//         uev = U(i,j+1,k);
//         uwv = U(i-1,j+1,k);
//         vn  = V(i,j+1,k);
//         vs  = V(i,j,k);
//         wfv = W(i,j+1,k);
//         wbv = W(i,j+1,k-1);
//     }
//     if (W(i,j,k) > 0)
//     {
//         uew = U(i,j,k);
//         uww = U(i-1,j,k);
//         vnw = V(i,j,k);
//         vsw = V(i,j-1,k);
//         wf  = W(i,j,k);
//         wb  = W(i,j,k-1);
//     }
//     else
//     {
//         uew = U(i,j+1,k);
//         uww = U(i-1,j,k+1);
//         vnw = V(i,j,k+1);
//         vsw = V(i,j-1,k+1);
//         wf  = W(i,j,k+1);
//         wb  = W(i,j,k);
//     }
//     un = ( U(i,j+1,k) + U(i,j,k)   ) / 2.0;
//     us = ( U(i,j,k)   + U(i,j-1,k) ) / 2.0;
//     uf = ( U(i,j,k+1) + U(i,j,k)   ) / 2.0;
//     ub = ( U(i,j,k)   + U(i,j,k-1) ) / 2.0;
//
//     ve = ( V(i+1,j,k) + V(i,j,k)   ) / 2.0;
//     vw = ( V(i,j,k)   + V(i-1,j,k) ) / 2.0;
//     vf = ( V(i,j,k+1) + V(i,j,k)   ) / 2.0;
//     vb = ( V(i,j,k)   + V(i,j,k-1) ) / 2.0;
//
//     we = ( W(i+1,j,k) + W(i,j,k)   ) / 2.0;
//     ww = ( W(i,j,k)   + W(i-1,j,k) ) / 2.0;
//     wn = ( W(i,j+1,k) + W(i,j,k)   ) / 2.0;
//     ws = ( W(i,j,k)   + W(i,j-1,k) ) / 2.0;
// }
