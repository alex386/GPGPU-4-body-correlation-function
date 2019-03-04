/*---------------------------------------*/
/* defines for float/double              */
/*---------------------------------------*/
#define BLOCK_SIZE 500


//single precision switch
#ifdef PSINGLE
#define real float
#define real4 float4
#define real3 float3
#define make_real3 make_float3
#define make_real4 make_float4
#endif

//double precision switch
#ifdef PDOUBLE
#define real double
#define real4 double4
#define real3 double3
#define make_real3 make_double3
#define make_real4 make_double4
#endif

