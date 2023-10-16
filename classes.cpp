#include "header_classes.hpp"

void param_diff(P_gas & gas)
{
  gas.segm_T = 1;
  gas.segm_X = 1;
  gas.mu = 0.1;
  gas.p_ro = 10;
  gas.p_gamma = 1;
}

void param_she (P_gas & gas, P_she & shem)
{
  shem.h_x = 0.01;
  shem.tau = 0.01;
  shem.m_x = (int) (gas.segm_X / shem.h_x);
  shem.n = (int) (gas.segm_T/shem.tau);
  shem.dim = shem.m_x + 1;
  shem.eta = gas.mu;
}

/*
        int parse_command_line (int argc, char * argv [])
        {
         if(!      ((argc==4)&&
             (sscanf(argv[1],"%le",&segm_T) == 1)&&
             (sscanf(argv[2],"%le",&segm_X) == 1)&&
             (sscanf(argv[3],"%le",&mu) == 1)&&
             (sscanf(argv[4],"%le",&x1) == 1)&&
             (sscanf(argv[5],"%le",&t0) == 1)&&
 		    (sscanf(argv[6],"%le",&t1) == 1)&&
 		    (sscanf(argv[7],"%d",&gamma_or_const) == 1)&&
 		    (sscanf(argv[8],"%le",&constant) == 1)&&
 		    (sscanf(argv[9],"%le",&viscosity) == 1)&&
             (gamma_or_const == 1))
                 ||((argc==8)&&
             (sscanf(argv[1],"%le",&m) == 1)&&
             (sscanf(argv[2],"%le",&n) == 1)&&
             (sscanf(argv[3],"%le",&x0) == 1)&&
             (sscanf(argv[4],"%le",&x1) == 1)&&
             (sscanf(argv[5],"%le",&t0) == 1)&&
 		    (sscanf(argv[6],"%le",&t1) == 1)&&
 		    (sscanf(argv[7],"%le",&viscosity) == 1))
                 || (argc < 8)
         )
       		{
       			return -1;
       		}
         return 0;
        }
*/
