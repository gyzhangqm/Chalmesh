#define R0 ((i_loc-1) * info->dr)
#define R1 (i_loc * info->dr)
#define ALPHA_0(r)   ((R1 - r)/info->dr)
#define D_ALPHA_0(r) (-1.0/info->dr)
#define ALPHA_1(r)   ((r - R0)/info->dr)
#define D_ALPHA_1(r) (1.0/info->dr)

#define S0 ((j_loc-1)  * info->ds)
#define S1 (j_loc * info->ds)
#define BETA_0(s)   ((S1 - s)/info->ds)
#define D_BETA_0(s) (-1.0/info->ds)
#define BETA_1(s)   ((s - S0)/info->ds)
#define D_BETA_1(s) (1.0/info->ds)

#define LIN(r,s,x) (ALPHA_0(r) * BETA_0(s) * x(i_loc  ,j_loc  ) +  \
		    ALPHA_1(r) * BETA_0(s) * x(i_loc+1,j_loc  ) +  \
		    ALPHA_0(r) * BETA_1(s) * x(i_loc  ,j_loc+1) +  \
		    ALPHA_1(r) * BETA_1(s) * x(i_loc+1,j_loc+1) )
#define LIN_R(r,s,x) (D_ALPHA_0(r) * BETA_0(s) * x(i_loc  ,j_loc  ) +  \
		      D_ALPHA_1(r) * BETA_0(s) * x(i_loc+1,j_loc  ) +  \
		      D_ALPHA_0(r) * BETA_1(s) * x(i_loc  ,j_loc+1) +  \
		      D_ALPHA_1(r) * BETA_1(s) * x(i_loc+1,j_loc+1) )
#define LIN_S(r,s,x) (ALPHA_0(r) * D_BETA_0(s) * x(i_loc  ,j_loc  ) +  \
		      ALPHA_1(r) * D_BETA_0(s) * x(i_loc+1,j_loc  ) +  \
		      ALPHA_0(r) * D_BETA_1(s) * x(i_loc  ,j_loc+1) +  \
		      ALPHA_1(r) * D_BETA_1(s) * x(i_loc+1,j_loc+1) )
#define LIN_RS(r,s,x) (D_ALPHA_0(r) * D_BETA_0(s) * x(i_loc  ,j_loc  ) +  \
		      D_ALPHA_1(r) * D_BETA_0(s) * x(i_loc+1,j_loc  ) +  \
		      D_ALPHA_0(r) * D_BETA_1(s) * x(i_loc  ,j_loc+1) +  \
		      D_ALPHA_1(r) * D_BETA_1(s) * x(i_loc+1,j_loc+1) )
