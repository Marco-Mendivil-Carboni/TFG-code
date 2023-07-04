//Structs

struct knot_type
{
  int invariant;
  char name[8];
};

struct crossing_info
{
  int hand;
  int i_o;
  int i_u;
  float t_o;
  float t_u;
  int i;
  int j;
  int k;
};

typedef struct crossing_info crossing_info;
typedef struct knot_type knot_type;

//Lookup table

const knot_type known_knot[] = 
{
  {     905463,    "3_1"},
  {    2509099,    "4_1"},
  {    2545745,    "5_1"},
  {    4925488,    "5_2"},
  {    8132760,    "6_1"},
  {   12240588,    "6_2"},
  {   17066075,    "6_3"},
  {    5080627,    "7_1"},
  {   12160074,    "7_2"},
  {   17161433,    "7_3"},
  {   22609223,    "7_4"},
  {   29272665,    "7_5"},
  {   36411894,    "7_6"},
  {   44444654,    "7_7"},
  {   16970983,    "8_1"},
  {   29649143,    "8_2"},
  {   29023769,    "8_3"},
  {   36551120,    "8_4"},
  {   45062872,    "8_5"},
  {   53487839,    "8_6"},
  {   53996320,    "8_7"},
  {   63138814,    "8_8"},
  {   63691163,    "8_9"},
  {   74235537,   "8_10"},
  {   73639120,   "8_11"},
  {   84681481,   "8_12"},
  {   84893732,   "8_13"},
  {   97004963,   "8_14"},
  {  110104951,   "8_15"},
  {  124614840,   "8_16"},
  {  139135171,   "8_17"},
  {  205711733,   "8_18"},
  {     972667,   "8_19"},
  {    8198629,   "8_20"},
  {   22718960,   "8_21"},
  {    8602498,    "9_1"},
  {   22609223,    "9_2"},
  {   37253791,    "9_3"},
  {   44752528,    "9_4"},
  {   53151207,    "9_5"},
  {   74834358,    "9_6"},
  {   85106248,    "9_7"},
  {   97004963,    "9_8"},
  {   98376073,    "9_9"},
  {  110346957,   "9_10"},
  {  111076775,   "9_11"},
  {  123585369,   "9_12"},
  {  138589445,   "9_13"},
  {  138047253,   "9_14"},
  {  153380336,   "9_15"},
  {  155390394,   "9_16"},
  {  154813639,   "9_17"},
  {  170046495,   "9_18"},
  {  169445857,   "9_19"},
  {  170952185,   "9_20"},
  {  186389865,   "9_21"},
  {  187969550,   "9_22"},
  {  204718107,   "9_23"},
  {  205711733,   "9_24"},
  {  222958011,   "9_25"},
  {  224340024,   "9_26"},
  {  243685842,   "9_27"},
  {  263925059,   "9_28"},
  {  263925059,   "9_29"},
  {  284874514,   "9_30"},
  {  306724656,   "9_31"},
  {  353171871,   "9_32"},
  {  377343177,   "9_33"},
  {  482280901,   "9_34"},
  {   73244041,   "9_35"},
  {  139407199,   "9_36"},
  {  204059022,   "9_37"},
  {  328437926,   "9_38"},
  {  305108326,   "9_39"},
  {  570041235,   "9_40"},
  {  242245383,   "9_41"},
  {    4976778,   "9_42"},
  {   17449975,   "9_43"},
  {   29148084,   "9_44"},
  {   53319390,   "9_45"},
  {    8132760,   "9_46"},
  {   74434274,   "9_47"},
  {   73441448,   "9_48"},
  {   63322107,   "9_49"},
  {   29023769,   "10_1"},
  {   55369240,   "10_2"},
  {   62773025,   "10_3"},
  {   73837058,   "10_4"},
  {  113042166,   "10_5"},
  {  140227359,   "10_6"},
  {  186704696,   "10_7"},
  {   86390821,   "10_8"},
  {  157132370,   "10_9"},
  {  204388432,  "10_10"},
  {  187019793,  "10_11"},
  {  225380122,  "10_12"},
  {  282929039,  "10_13"},
  {  331376270,  "10_14"},
  {  188921713,  "10_15"},
  {  223302331,  "10_16"},
  {  173388341,  "10_17"},
  {  305511092,  "10_18"},
  {  265428501,  "10_19"},
  {  123841753,  "10_20"},
  {  207039303,  "10_21"},
  {  244769809,  "10_22"},
  {  354476575,  "10_23"},
  {  305511092,  "10_24"},
  {  430119188,  "10_25"},
  {  378691749,  "10_26"},
  {  512679691,  "10_27"},
  {  283705017,  "10_28"},
  {  402892359,  "10_29"},
  {  452933450,  "10_30"},
  {  328020316,  "10_31"},
  {  484312489,  "10_32"},
  {  426294598,  "10_33"},
  {  138318216,  "10_34"},
  {  241886750,  "10_35"},
  {  262425887,  "10_36"},
  {  283705017,  "10_37"},
  {  351437316,  "10_38"},
  {  379140448,  "10_39"},
  {  571698481,  "10_40"},
  {  511110381,  "10_41"},
  {  664390717,  "10_42"},
  {  540107174,  "10_43"},
  {  632186650,  "10_44"},
  {  801532507,  "10_45"},
  {   99763102,  "10_46"},
  {  173691998,  "10_47"},
  {  246954867,  "10_48"},
  {  355783685,  "10_49"},
  {  286436402,  "10_50"},
  {  456382897,  "10_51"},
  {  354476575,  "10_52"},
  {  538496402,  "10_53"},
  {  225380122,  "10_54"},
  {  375997010,  "10_55"},
  {  430119188,  "10_56"},
  {  633931832,  "10_57"},
  {  425818806,  "10_58"},
  {  570041235,  "10_59"},
  {  731354331,  "10_60"},
  {  111565395,  "10_61"},
  {  208716151,  "10_62"},
  {  328437926,  "10_63"},
  {  267326687,  "10_64"},
  {  403822455,  "10_65"},
  {  573358132,  "10_66"},
  {  400578102,  "10_67"},
  {  328020316,  "10_68"},
  {  766121168,  "10_69"},
  {  455394089,  "10_70"},
  {  600641665,  "10_71"},
  {  542256987,  "10_72"},
  {  697546628,  "10_73"},
  {  400578102,  "10_74"},
  {  664390717,  "10_75"},
  {  330956795,  "10_76"},
  {  403822455,  "10_77"},
  {  482787246,  "10_78"},
  {  381408425,  "10_79"},
  {  514251407,  "10_80"},
  {  731977836,  "10_81"},
  {  408027545,  "10_82"},
  {  699989471,  "10_83"},
  {  768681169,  "10_84"},
  {  334764675,  "10_85"},
  {  733855620,  "10_86"},
  {  666774847,  "10_87"},
  { 1032094805,  "10_88"},
  {  991857978,  "10_89"},
  {  602342779,  "10_90"},
  {  546589063,  "10_91"},
  {  804804748,  "10_92"},
  {  456875464,  "10_93"},
  {  517416399,  "10_94"},
  {  840587430,  "10_95"},
  {  874925245,  "10_96"},
  {  763565437,  "10_97"},
  {  666774847,  "10_98"},
  {  672175175,  "10_99"},
  {  434458686, "10_100"},
  {  730102457, "10_101"},
  {  541720353, "10_102"},
  {  571698481, "10_103"},
  {  607476095, "10_104"},
  {  838577642, "10_105"},
  {  576699794, "10_106"},
  {  875607196, "10_107"},
  {  404285798, "10_108"},
  {  739520554, "10_109"},
  {  698155554, "10_110"},
  {  602908636, "10_111"},
  {  776410202, "10_112"},
  { 1249966536, "10_113"},
  {  878343868, "10_114"},
  { 1202239772, "10_115"},
  {  924843234, "10_116"},
  { 1076349691, "10_117"},
  {  963710742, "10_118"},
  { 1035065795, "10_119"},
  { 1113614483, "10_120"},
  { 1341218615, "10_121"},
  { 1119019553, "10_122"},
  { 1498319904, "10_123"},
  {     148922, "10_124"},
  {   12484463, "10_125"},
  {   36831640, "10_126"},
  {   85747332, "10_127"},
  {   12730743, "10_128"},
  {   63138814, "10_129"},
  {   29272665, "10_130"},
  {   97004963, "10_131"},
  {    2545745, "10_132"},
  {   36411894, "10_133"},
  {   54507206, "10_134"},
  {  138318216, "10_135"},
  {   22718960, "10_136"},
  {   62955787, "10_137"},
  {  124872290, "10_138"},
  {     777715, "10_139"},
  {    8198629, "10_140"},
  {   45062872, "10_141"},
  {   23384986, "10_142"},
  {   74235537, "10_143"},
  {  153665944, "10_144"},
  {     883662, "10_145"},
  {  109863211, "10_146"},
  {   73639120, "10_147"},
  {   97689315, "10_148"},
  {  170952185, "10_149"},
  {   85747332, "10_150"},
  {  187653654, "10_151"},
  {   11757893, "10_152"},
  {      79269, "10_153"},
  {   16780738, "10_154"},
  {   63691163, "10_155"},
  {  124614840, "10_156"},
  {  244045806, "10_157"},
  {  205381258, "10_158"},
  {  154526966, "10_159"},
  {   45217741, "10_160"},
  {    2436278, "10_161"},
  {  123841753, "10_162"},
  {  263925059, "10_163"},
  {  204388432, "10_164"},
  {  153380336, "10_165"}
};

//Functions

double calc_det( int N, double **M)
{
  int *P = (int*)malloc((N+1)*sizeof(int));
  for( int i = 0; i<(N+1); i++)
  {
    P[i] = i;
  }
  for( int j = 0; j<N; j++)
  {
    int i_max;
    double M_max = 0.0;
    for( int i = j; i<N; i++)
    {
      if( fabs(M[i][j])>M_max)
      { 
        i_max = i;
        M_max = fabs(M[i][j]);
      }
    }
    if( M_max<__FLT_MIN__){ return 0.0;}
    if( i_max!=j)
    {
      int P_aux;
      double *M_aux;
      P_aux = P[j];
      P[j] = P[i_max];
      P[i_max] = P_aux;
      M_aux = M[j];
      M[j] = M[i_max];
      M[i_max] = M_aux;
      P[N]++;
    }
    for( int i=j+1; i<N; i++)
    {
      M[i][j] /= M[j][j];
      for( int k = j+1; k<N; k++)
      {
        M[i][k] -= M[i][j]*M[j][k];
      }
    }
  }
  double det = 1.0;
  for( int i = 0; i<N; i++)
  {
    det *= M[i][i];
  }
  if( (P[N]-N)%2!=0)
  {
    det *= (-1.0);
  }
  return det;
}

void calc_invariant( int N, float *r, char *knot_name)
{
  //close chain
  float r_cm[3];
  r_cm[0] = r_cm[1] = r_cm[2] = 0.0;
  for( int i_p = 0; i_p<N; i_p++)
  {
    r_cm[0] += r[3*i_p+0];
    r_cm[1] += r[3*i_p+1];
    r_cm[2] += r[3*i_p+2];
  }
  r_cm[0] /= N;
  r_cm[1] /= N;
  r_cm[2] /= N;
  float dist_max = 0.0;
  for( int i_p = 0; i_p<N; i_p++)
  {
    float dist = 0.0;
    for( int i_c = 0; i_c<3; i_c++)
    {
      dist += (r[3*i_p+i_c]-r_cm[i_c])*(r[3*i_p+i_c]-r_cm[i_c]);
    }
    dist = sqrt(dist);
    if( dist>dist_max)
    {
      dist_max = dist;
    }
  }
  float *r_l = (float*)malloc(3*(N+4)*sizeof(float));
  for( int i_p = 0; i_p<N; i_p++)
  {
    for( int i_c = 0; i_c<3; i_c++)
    {
      r_l[3*i_p+i_c] = r[3*i_p+i_c];
    }
  }
  float norm_aux[3] = {0.0,0.0,0.0};
  for( int i_c = 0; i_c<3; i_c++)
  {
    r_l[3*(N)+i_c] = r[3*(N-1)+i_c]-r_cm[i_c];
    r_l[3*(N+2)+i_c] = r[3*(0)+i_c]-r_cm[i_c];
    norm_aux[0] += r_l[3*(N+0)+i_c]*r_l[3*(N+0)+i_c];
    norm_aux[2] += r_l[3*(N+2)+i_c]*r_l[3*(N+2)+i_c];
  }
  norm_aux[0] = sqrt(norm_aux[0]);
  norm_aux[2] = sqrt(norm_aux[2]);
  for( int i_c = 0; i_c<3; i_c++)
  {
    r_l[3*(N+0)+i_c] = 1.5*dist_max*r_l[3*(N+0)+i_c]/norm_aux[0];
    r_l[3*(N+2)+i_c] = 1.5*dist_max*r_l[3*(N+2)+i_c]/norm_aux[2];
    r_l[3*(N+1)+i_c] = r_l[3*(N+0)+i_c]+r_l[3*(N+2)+i_c];
    norm_aux[1] += r_l[3*(N+1)+i_c]*r_l[3*(N+1)+i_c];
  }
  norm_aux[1] = sqrt(norm_aux[1]);
  for( int i_c = 0; i_c<3; i_c++)
  {
    r_l[3*(N+1)+i_c] = 1.5*dist_max*r_l[3*(N+1)+i_c]/norm_aux[1];
    r_l[3*(N+0)+i_c] += r_cm[i_c];
    r_l[3*(N+1)+i_c] += r_cm[i_c];
    r_l[3*(N+2)+i_c] += r_cm[i_c];
    r_l[3*(N+3)+i_c] = r[3*(0)+i_c];
  }
  //find crossings
  int n_crossings = 0;
  crossing_info *crossing = NULL;
  for( int i_l = 0; i_l<(N+3); i_l++)
  {
    for( int j_l = i_l+2; j_l<(N+3); j_l++)
    {
      float link_i[3], link_j[3], sep_ij[3];
      for( int i_c = 0; i_c<3; i_c++)
      {
        link_i[i_c] = r_l[3*(i_l+1)+i_c]-r_l[3*i_l+i_c];
        link_j[i_c] = r_l[3*(j_l+1)+i_c]-r_l[3*j_l+i_c];
        sep_ij[i_c] = r_l[3*j_l+i_c]-r_l[3*i_l+i_c];
      }
      float cross_prod = link_i[0]*link_j[1]-link_i[1]*link_j[0];
      float t_i = (sep_ij[0]*link_j[1]-sep_ij[1]*link_j[0])/cross_prod;
      float t_j = (sep_ij[0]*link_i[1]-sep_ij[1]*link_i[0])/cross_prod;
      if( (0.0<t_i)&&(t_i<1.0)&&(0.0<t_j)&&(t_j<1.0))
      {
        n_crossings++;
        crossing = (crossing_info*)realloc(crossing,n_crossings*sizeof(crossing_info));
        if( (r_l[3*i_l+2]+t_i*link_i[2])>(r_l[3*j_l+2]+t_j*link_j[2]))
        {
          crossing[n_crossings-1].hand = (+1)*((cross_prod>0)-(cross_prod<0));
          crossing[n_crossings-1].i_o = i_l;
          crossing[n_crossings-1].i_u = j_l;
          crossing[n_crossings-1].t_o = t_i;
          crossing[n_crossings-1].t_u = t_j;
        }
        else
        {
          crossing[n_crossings-1].hand = (-1)*((cross_prod>0)-(cross_prod<0));
          crossing[n_crossings-1].i_o = j_l;
          crossing[n_crossings-1].i_u = i_l;
          crossing[n_crossings-1].t_o = t_j;
          crossing[n_crossings-1].t_u = t_i;
        }
      }
    }
  }
  free(r_l);
  //locate last underpass
  int idx_last_u;
  float t_last_u = -1.0;
  for( int i_l = N+2; i_l>0; i_l--)
  {
    for( int l = 0; l<n_crossings; l++)
    {
      if( (i_l==crossing[l].i_u)&&(crossing[l].t_u>t_last_u))
      {
        t_last_u = crossing[l].t_u;
        idx_last_u = l;
      }
    }
    if( t_last_u>0.0)
    {
      break;
    }
  }
  //number arcs
  int arc_idx = 0;
  int *idx_l = (int*)malloc(n_crossings*sizeof(int));
  float *t_l = (float*)malloc(n_crossings*sizeof(int));
  int *o_l = (int*)malloc(n_crossings*sizeof(int));
  for( int i_l = 0; i_l<(N+3); i_l++)
  {
    int n_l = 0;
    for( int l = 0; l<n_crossings; l++)
    {
      if( (crossing[l].i_o==i_l))
      {
        int j_l_new = n_l;
        for( int j_l = 0; j_l<n_l; j_l++)
        {
          if( (crossing[l].t_o<t_l[j_l]))
          {
            j_l_new = j_l; break;
          }
        }
        for( int j_l = n_l; j_l>j_l_new; j_l--)
        {
          idx_l[j_l] = idx_l[j_l-1];
          t_l[j_l] = t_l[j_l-1];
          o_l[j_l] = o_l[j_l-1];
        }
        idx_l[j_l_new] = l;
        t_l[j_l_new] = crossing[l].t_o;
        o_l[j_l_new] = 1;
        n_l++;
      }
      if( (crossing[l].i_u==i_l))
      {
        int j_l_new = n_l;
        for( int j_l = 0; j_l<n_l; j_l++)
        {
          if( (crossing[l].t_u<t_l[j_l]))
          {
            j_l_new = j_l; break;
          }
        }
        for( int j_l = n_l; j_l>j_l_new; j_l--)
        {
          idx_l[j_l] = idx_l[j_l-1];
          t_l[j_l] = t_l[j_l-1];
          o_l[j_l] = o_l[j_l-1];
        }
        idx_l[j_l_new] = l;
        t_l[j_l_new] = crossing[l].t_u;
        o_l[j_l_new] = 0;
        n_l++;
      }
    }
    for( int j_l = 0; j_l<n_l; j_l++)
    {
      if( o_l[j_l])
      {
        crossing[idx_l[j_l]].i = arc_idx;
      }
      else
      {
        crossing[idx_l[j_l]].j = arc_idx;
        if( idx_l[j_l]!=idx_last_u){ arc_idx++;}else{ arc_idx=0;}
        crossing[idx_l[j_l]].k = arc_idx;
      }
    }
  }
  free(idx_l);
  free(t_l);
  free(o_l);
  //allocate matrices
  int ***A = (int***)malloc(n_crossings*sizeof(int**));
  for( int i = 0; i<n_crossings; i++)
  {
    A[i] = (int**)malloc(n_crossings*sizeof(int*));
    for( int j = 0; j<n_crossings; j++)
    {
      A[i][j] = (int*)calloc(2,sizeof(int));
    }
  }
  int n_M = n_crossings-1;
  double **M = (double**)malloc(n_M*sizeof(double*));
  for( int i = 0; i<n_M; i++)
  {
    M[i] = (double*)malloc(n_M*sizeof(double));
  }
  //make Alexander matrix
  for( int l = 0; l<n_crossings; l++)
  {
    if( crossing[l].hand==1)
    {
      A[l][crossing[l].i][0] += (+1);
      A[l][crossing[l].i][1] += (-1);
      A[l][crossing[l].j][0] += (-1);
      A[l][crossing[l].k][1] += (+1);
    }
    else
    {
      A[l][crossing[l].i][0] += (+1);
      A[l][crossing[l].i][1] += (-1);
      A[l][crossing[l].j][1] += (+1);
      A[l][crossing[l].k][0] += (-1);
    }
  }
  //calculate knot invariant
  double fp_k_inv = 1.0;
  double t = -1.1;
  for( int i = 0; i<n_M; i++)
  {
    for( int j = 0; j<n_M; j++)
    {
      M[i][j] = A[i][j][0]+t*A[i][j][1];
    }
  }
  fp_k_inv *= calc_det(n_M,M);
  for( int i = 0; i<n_M; i++)
  {
    for( int j = 0; j<n_M; j++)
    {
      M[i][j] = A[i][j][0]+(1.0/t)*A[i][j][1];
    }
  }
  fp_k_inv *= calc_det(n_M,M);
  //deallocate matrices
  for( int i = 0; i<n_M; i++)
  {
    free(M[i]);
  }
  free(M);
  for( int i = 0; i<n_crossings; i++)
  {
    for( int j = 0; j<n_crossings; j++)
    {
      free(A[i][j]);
    }
    free(A[i]);
  }
  free(A);
  //return knot name
  int int_k_inv = (int)round(100000.0*fp_k_inv);
  if( int_k_inv==0)
  {
    snprintf(knot_name,sizeof(knot_name),"f");
    return;
  }
  if( int_k_inv==100000)
  {
    snprintf(knot_name,sizeof(knot_name),"0");
    return;
  }
  for( int i = 0; i<249; i++)
  {
    if( int_k_inv==known_knot[i].invariant)
    {
      snprintf(knot_name,sizeof(knot_name),"%s",known_knot[i].name);
      return;
    }
  }
  snprintf(knot_name,sizeof(knot_name),"%d",int_k_inv);
}
