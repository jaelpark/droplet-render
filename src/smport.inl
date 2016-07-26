#ifndef SMPORT_INL
#define SMPORT_INL

//Some stuff ported from DirectXMath - temporary placeholders for some other standard-compatible library.
//Probably have to dump these anyway, at least when the AVX implementation comes.
//------------------------------------------------------------------------------

const float XM_PI           = 3.141592654f;
const float XM_2PI          = 6.283185307f;
const float XM_1DIVPI       = 0.318309886f;
const float XM_1DIV2PI      = 0.159154943f;
const float XM_PIDIV2       = 1.570796327f;
const float XM_PIDIV4       = 0.785398163f;

#define XMVECTOR __m128
#define FXMVECTOR __m128
#define XM_PERMUTE_PS(v, c) _mm_shuffle_ps(v,v,c)

#define g_XMMaskY _mm_castsi128_ps(_mm_set_epi32(0,0,-1,0))
#define g_XMMask3 _mm_castsi128_ps(_mm_set_epi32(0,-1,-1,-1))

typedef struct XMMATRIX{
    XMVECTOR r[4];
} FXMMATRIX __attribute__((aligned(16)));

//------------------------------------------------------------------------------

inline XMVECTOR XMVectorMultiplyAdd
(
    FXMVECTOR V1,
    FXMVECTOR V2,
    FXMVECTOR V3
)
{
    XMVECTOR vResult = _mm_mul_ps( V1, V2 );
    return _mm_add_ps(vResult, V3 );
}

inline XMVECTOR XMVector3Dot
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(USE_SSE4)
    return _mm_dp_ps(V1,V2,0x7f);
#else
    /*// Perform the dot product
    XMVECTOR vDot = _mm_mul_ps(V1,V2);
    // x=Dot.vector4_f32[1], y=Dot.vector4_f32[2]
    XMVECTOR vTemp = XM_PERMUTE_PS(vDot,_MM_SHUFFLE(2,1,2,1));
    // Result.vector4_f32[0] = x+y
    vDot = _mm_add_ss(vDot,vTemp);
    // x=Dot.vector4_f32[2]
    vTemp = XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(1,1,1,1));
    // Result.vector4_f32[0] = (x+y)+z
    vDot = _mm_add_ss(vDot,vTemp);
    // Splat x
    return XM_PERMUTE_PS(vDot,_MM_SHUFFLE(0,0,0,0));*/
    XMVECTOR vTemp = _mm_mul_ps(V1,V2);
    vTemp = _mm_and_ps( vTemp, g_XMMask3 );
    vTemp = _mm_hadd_ps(vTemp,vTemp);
    return _mm_hadd_ps(vTemp,vTemp);
#endif
}

inline XMVECTOR XMVector3Length
(
    FXMVECTOR V
)
{
    XMVECTOR d = XMVector3Dot(V,V);
    return _mm_sqrt_ps(d);
}

inline XMVECTOR XMVector3Normalize
(
    FXMVECTOR V
)
{
    XMVECTOR l = XMVector3Length(V);
    return _mm_div_ps(V,l);
}

inline XMVECTOR XMVector4Dot
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
#if defined(USE_SSE4)
    return _mm_dp_ps(V1,V2,0xff);
#else
    /*XMVECTOR vTemp2 = V2;
    XMVECTOR vTemp = _mm_mul_ps(V1,vTemp2);
    vTemp2 = _mm_shuffle_ps(vTemp2,vTemp,_MM_SHUFFLE(1,0,0,0)); // Copy X to the Z position and Y to the W position
    vTemp2 = _mm_add_ps(vTemp2,vTemp);          // Add Z = X+Z; W = Y+W;
    vTemp = _mm_shuffle_ps(vTemp,vTemp2,_MM_SHUFFLE(0,3,0,0));  // Copy W to the Z position
    vTemp = _mm_add_ps(vTemp,vTemp2);           // Add Z and W together
    return XM_PERMUTE_PS(vTemp,_MM_SHUFFLE(2,2,2,2));    // Splat Z and return*/
    XMVECTOR r = _mm_mul_ps(V1,V2);
    r = _mm_hadd_ps(r,r);
    return _mm_hadd_ps(r,r);
#endif
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3Cross
(
    FXMVECTOR V1,
    FXMVECTOR V2
)
{
    // [ V1.y*V2.z - V1.z*V2.y, V1.z*V2.x - V1.x*V2.z, V1.x*V2.y - V1.y*V2.x ]
    // y1,z1,x1,w1
    XMVECTOR vTemp1 = XM_PERMUTE_PS(V1,_MM_SHUFFLE(3,0,2,1));
    // z2,x2,y2,w2
    XMVECTOR vTemp2 = XM_PERMUTE_PS(V2,_MM_SHUFFLE(3,1,0,2));
    // Perform the left operation
    XMVECTOR vResult = _mm_mul_ps(vTemp1,vTemp2);
    // z1,x1,y1,w1
    vTemp1 = XM_PERMUTE_PS(vTemp1,_MM_SHUFFLE(3,0,2,1));
    // y2,z2,x2,w2
    vTemp2 = XM_PERMUTE_PS(vTemp2,_MM_SHUFFLE(3,1,0,2));
    // Perform the right operation
    vTemp1 = _mm_mul_ps(vTemp1,vTemp2);
    // Subract the right from left, and return answer
    vResult = _mm_sub_ps(vResult,vTemp1);
    // Set w to zero
    return _mm_and_ps(vResult,g_XMMask3);
}

inline XMVECTOR XMVectorNegate
(
    FXMVECTOR V
)
{
    XMVECTOR Z;
    Z = _mm_setzero_ps();
    return _mm_sub_ps( Z, V );
}

//------------------------------------------------------------------------------

inline XMVECTOR XMVector3TransformCoord
(
    FXMVECTOR V,
    FXMMATRIX M
)
{
    XMVECTOR Z = FL_PERMUTE(V,_MM_SHUFFLE(2,2,2,2));//XMVectorSplatZ(V);
    XMVECTOR Y = FL_PERMUTE(V,_MM_SHUFFLE(1,1,1,1));//XMVectorSplatY(V);
    XMVECTOR X = FL_PERMUTE(V,_MM_SHUFFLE(0,0,0,0));//XMVectorSplatX(V);

    XMVECTOR Result = XMVectorMultiplyAdd(Z, M.r[2], M.r[3]);
    Result = XMVectorMultiplyAdd(Y, M.r[1], Result);
    Result = XMVectorMultiplyAdd(X, M.r[0], Result);

    XMVECTOR W = FL_PERMUTE(Result,_MM_SHUFFLE(3,3,3,3));//XMVectorSplatW(Result);
    return _mm_div_ps(Result,W);
    //return XMVectorDivide( Result, W );
}

//------------------------------------------------------------------------------

inline XMVECTOR XMQuaternionConjugate
(
    FXMVECTOR Q
)
{
	//return _mm_mul_ps(Q,_mm_set_ps(-1.0f,-1.0f,-1.0f,1.0f));
	return _mm_mul_ps(Q,float4(-1.0f,-1.0f,-1.0f,1.0f).v);
}

inline XMVECTOR XMQuaternionMultiply
(
    FXMVECTOR Q1,
    FXMVECTOR Q2
)
{
    // Returns the product Q2*Q1 (which is the concatenation of a rotation Q1 followed by the rotation Q2)

    // [ (Q2.w * Q1.x) + (Q2.x * Q1.w) + (Q2.y * Q1.z) - (Q2.z * Q1.y),
    //   (Q2.w * Q1.y) - (Q2.x * Q1.z) + (Q2.y * Q1.w) + (Q2.z * Q1.x),
    //   (Q2.w * Q1.z) + (Q2.x * Q1.y) - (Q2.y * Q1.x) + (Q2.z * Q1.w),
    //   (Q2.w * Q1.w) - (Q2.x * Q1.x) - (Q2.y * Q1.y) - (Q2.z * Q1.z) ]

    /*static const XMVECTORF32 ControlWZYX = { 1.0f,-1.0f, 1.0f,-1.0f};
    static const XMVECTORF32 ControlZWXY = { 1.0f, 1.0f,-1.0f,-1.0f};
    static const XMVECTORF32 ControlYXWZ = {-1.0f, 1.0f, 1.0f,-1.0f};*/
    // Copy to SSE registers and use as few as possible for x86
    XMVECTOR Q2X = Q2;
    XMVECTOR Q2Y = Q2;
    XMVECTOR Q2Z = Q2;
    XMVECTOR vResult = Q2;
    // Splat with one instruction
    vResult = XM_PERMUTE_PS(vResult,_MM_SHUFFLE(3,3,3,3));
    Q2X = XM_PERMUTE_PS(Q2X,_MM_SHUFFLE(0,0,0,0));
    Q2Y = XM_PERMUTE_PS(Q2Y,_MM_SHUFFLE(1,1,1,1));
    Q2Z = XM_PERMUTE_PS(Q2Z,_MM_SHUFFLE(2,2,2,2));
    // Retire Q1 and perform Q1*Q2W
    vResult = _mm_mul_ps(vResult,Q1);
    XMVECTOR Q1Shuffle = Q1;
    // Shuffle the copies of Q1
    Q1Shuffle = XM_PERMUTE_PS(Q1Shuffle,_MM_SHUFFLE(0,1,2,3));
    // Mul by Q1WZYX
    Q2X = _mm_mul_ps(Q2X,Q1Shuffle);
    Q1Shuffle = XM_PERMUTE_PS(Q1Shuffle,_MM_SHUFFLE(2,3,0,1));
    // Flip the signs on y and z
	//Q2X = _mm_mul_ps(Q2X,_mm_set_ps(1.0f,-1.0f, -1.0f,1.0f));
	Q2X = _mm_mul_ps(Q2X,float4(1.0f,-1.0f, 1.0f,-1.0f).v);
    // Mul by Q1ZWXY
    Q2Y = _mm_mul_ps(Q2Y,Q1Shuffle);
    Q1Shuffle = XM_PERMUTE_PS(Q1Shuffle,_MM_SHUFFLE(0,1,2,3));
    // Flip the signs on z and w
	//Q2Y = _mm_mul_ps(Q2Y,_mm_set_ps(1.0f, 1.0f,-1.0f,-1.0f));
	Q2Y = _mm_mul_ps(Q2Y,float4(1.0f, 1.0f,-1.0f,-1.0f).v);
    // Mul by Q1YXWZ
    Q2Z = _mm_mul_ps(Q2Z,Q1Shuffle);
    vResult = _mm_add_ps(vResult,Q2X);
    // Flip the signs on x and w
	//Q2Z = _mm_mul_ps(Q2Z,_mm_set_ps(-1.0f, 1.0f, 1.0f,-1.0f));
	Q2Z = _mm_mul_ps(Q2Z,float4(-1.0f, 1.0f, 1.0f,-1.0f).v);
    Q2Y = _mm_add_ps(Q2Y,Q2Z);
    vResult = _mm_add_ps(vResult,Q2Y);
    return vResult;
}

//------------------------------------------------------------------------------
// Transform a vector using a rotation expressed as a unit quaternion

inline XMVECTOR XMVector3Rotate
(
    FXMVECTOR V,
    FXMVECTOR RotationQuaternion
)
{
    //XMVECTOR A = XMVectorSelect(g_XMSelect1110.v, V, g_XMSelect1110.v);
    XMVECTOR A = float4::select(float4::zero(),float4(V),float4::selectctrl(1,1,1,0)).v;

    XMVECTOR Q = XMQuaternionConjugate(RotationQuaternion);
    XMVECTOR Result = XMQuaternionMultiply(Q, A);
    return XMQuaternionMultiply(Result, RotationQuaternion);
}

//------------------------------------------------------------------------------

inline void XMScalarSinCos
(
    float* pSin,
    float* pCos,
    float  Value
)
{

    // Map Value to y in [-pi,pi], x = 2*pi*quotient + remainder.
    float quotient = XM_1DIV2PI*Value;
    if (Value >= 0.0f)
    {
        quotient = (float)((int)(quotient + 0.5f));
    }
    else
    {
        quotient = (float)((int)(quotient - 0.5f));
    }
    float y = Value - XM_2PI*quotient;

    // Map y to [-pi/2,pi/2] with sin(y) = sin(Value).
    float sign;
    if (y > XM_PIDIV2)
    {
        y = XM_PI - y;
        sign = -1.0f;
    }
    else if (y < -XM_PIDIV2)
    {
        y = -XM_PI - y;
        sign = -1.0f;
    }
    else
    {
        sign = +1.0f;
    }

    float y2 = y * y;

    // 11-degree minimax approximation
    *pSin = ( ( ( ( (-2.3889859e-08f * y2 + 2.7525562e-06f) * y2 - 0.00019840874f ) * y2 + 0.0083333310f ) * y2 - 0.16666667f ) * y2 + 1.0f ) * y;

    // 10-degree minimax approximation
    float p = ( ( ( ( -2.6051615e-07f * y2 + 2.4760495e-05f ) * y2 - 0.0013888378f ) * y2 + 0.041666638f ) * y2 - 0.5f ) * y2 + 1.0f;
    *pCos = sign*p;
}

inline XMMATRIX XMMatrixTranspose
(
    FXMMATRIX M
)
{
    // x.x,x.y,y.x,y.y
    XMVECTOR vTemp1 = _mm_shuffle_ps(M.r[0],M.r[1],_MM_SHUFFLE(1,0,1,0));
    // x.z,x.w,y.z,y.w
    XMVECTOR vTemp3 = _mm_shuffle_ps(M.r[0],M.r[1],_MM_SHUFFLE(3,2,3,2));
    // z.x,z.y,w.x,w.y
    XMVECTOR vTemp2 = _mm_shuffle_ps(M.r[2],M.r[3],_MM_SHUFFLE(1,0,1,0));
    // z.z,z.w,w.z,w.w
    XMVECTOR vTemp4 = _mm_shuffle_ps(M.r[2],M.r[3],_MM_SHUFFLE(3,2,3,2));
    XMMATRIX mResult;

    // x.x,y.x,z.x,w.x
    mResult.r[0] = _mm_shuffle_ps(vTemp1, vTemp2,_MM_SHUFFLE(2,0,2,0));
    // x.y,y.y,z.y,w.y
    mResult.r[1] = _mm_shuffle_ps(vTemp1, vTemp2,_MM_SHUFFLE(3,1,3,1));
    // x.z,y.z,z.z,w.z
    mResult.r[2] = _mm_shuffle_ps(vTemp3, vTemp4,_MM_SHUFFLE(2,0,2,0));
    // x.w,y.w,z.w,w.w
    mResult.r[3] = _mm_shuffle_ps(vTemp3, vTemp4,_MM_SHUFFLE(3,1,3,1));
    return mResult;
}

//------------------------------------------------------------------------------
// Return the inverse and the determinant of a 4x4 matrix
inline XMMATRIX XMMatrixInverse
(
    XMVECTOR* pDeterminant,
    FXMMATRIX  M
)
{
    XMMATRIX MT = XMMatrixTranspose(M);
    XMVECTOR V00 = XM_PERMUTE_PS(MT.r[2],_MM_SHUFFLE(1,1,0,0));
    XMVECTOR V10 = XM_PERMUTE_PS(MT.r[3],_MM_SHUFFLE(3,2,3,2));
    XMVECTOR V01 = XM_PERMUTE_PS(MT.r[0],_MM_SHUFFLE(1,1,0,0));
    XMVECTOR V11 = XM_PERMUTE_PS(MT.r[1],_MM_SHUFFLE(3,2,3,2));
    XMVECTOR V02 = _mm_shuffle_ps(MT.r[2], MT.r[0],_MM_SHUFFLE(2,0,2,0));
    XMVECTOR V12 = _mm_shuffle_ps(MT.r[3], MT.r[1],_MM_SHUFFLE(3,1,3,1));

    XMVECTOR D0 = _mm_mul_ps(V00,V10);
    XMVECTOR D1 = _mm_mul_ps(V01,V11);
    XMVECTOR D2 = _mm_mul_ps(V02,V12);

    V00 = XM_PERMUTE_PS(MT.r[2],_MM_SHUFFLE(3,2,3,2));
    V10 = XM_PERMUTE_PS(MT.r[3],_MM_SHUFFLE(1,1,0,0));
    V01 = XM_PERMUTE_PS(MT.r[0],_MM_SHUFFLE(3,2,3,2));
    V11 = XM_PERMUTE_PS(MT.r[1],_MM_SHUFFLE(1,1,0,0));
    V02 = _mm_shuffle_ps(MT.r[2],MT.r[0],_MM_SHUFFLE(3,1,3,1));
    V12 = _mm_shuffle_ps(MT.r[3],MT.r[1],_MM_SHUFFLE(2,0,2,0));

    V00 = _mm_mul_ps(V00,V10);
    V01 = _mm_mul_ps(V01,V11);
    V02 = _mm_mul_ps(V02,V12);
    D0 = _mm_sub_ps(D0,V00);
    D1 = _mm_sub_ps(D1,V01);
    D2 = _mm_sub_ps(D2,V02);
    // V11 = D0Y,D0W,D2Y,D2Y
    V11 = _mm_shuffle_ps(D0,D2,_MM_SHUFFLE(1,1,3,1));
    V00 = XM_PERMUTE_PS(MT.r[1], _MM_SHUFFLE(1,0,2,1));
    V10 = _mm_shuffle_ps(V11,D0,_MM_SHUFFLE(0,3,0,2));
    V01 = XM_PERMUTE_PS(MT.r[0], _MM_SHUFFLE(0,1,0,2));
    V11 = _mm_shuffle_ps(V11,D0,_MM_SHUFFLE(2,1,2,1));
    // V13 = D1Y,D1W,D2W,D2W
    XMVECTOR V13 = _mm_shuffle_ps(D1,D2,_MM_SHUFFLE(3,3,3,1));
    V02 = XM_PERMUTE_PS(MT.r[3], _MM_SHUFFLE(1,0,2,1));
    V12 = _mm_shuffle_ps(V13,D1,_MM_SHUFFLE(0,3,0,2));
    XMVECTOR V03 = XM_PERMUTE_PS(MT.r[2],_MM_SHUFFLE(0,1,0,2));
    V13 = _mm_shuffle_ps(V13,D1,_MM_SHUFFLE(2,1,2,1));

    XMVECTOR C0 = _mm_mul_ps(V00,V10);
    XMVECTOR C2 = _mm_mul_ps(V01,V11);
    XMVECTOR C4 = _mm_mul_ps(V02,V12);
    XMVECTOR C6 = _mm_mul_ps(V03,V13);

    // V11 = D0X,D0Y,D2X,D2X
    V11 = _mm_shuffle_ps(D0,D2,_MM_SHUFFLE(0,0,1,0));
    V00 = XM_PERMUTE_PS(MT.r[1], _MM_SHUFFLE(2,1,3,2));
    V10 = _mm_shuffle_ps(D0,V11,_MM_SHUFFLE(2,1,0,3));
    V01 = XM_PERMUTE_PS(MT.r[0], _MM_SHUFFLE(1,3,2,3));
    V11 = _mm_shuffle_ps(D0,V11,_MM_SHUFFLE(0,2,1,2));
    // V13 = D1X,D1Y,D2Z,D2Z
    V13 = _mm_shuffle_ps(D1,D2,_MM_SHUFFLE(2,2,1,0));
    V02 = XM_PERMUTE_PS(MT.r[3], _MM_SHUFFLE(2,1,3,2));
    V12 = _mm_shuffle_ps(D1,V13,_MM_SHUFFLE(2,1,0,3));
    V03 = XM_PERMUTE_PS(MT.r[2],_MM_SHUFFLE(1,3,2,3));
    V13 = _mm_shuffle_ps(D1,V13,_MM_SHUFFLE(0,2,1,2));

    V00 = _mm_mul_ps(V00,V10);
    V01 = _mm_mul_ps(V01,V11);
    V02 = _mm_mul_ps(V02,V12);
    V03 = _mm_mul_ps(V03,V13);
    C0 = _mm_sub_ps(C0,V00);
    C2 = _mm_sub_ps(C2,V01);
    C4 = _mm_sub_ps(C4,V02);
    C6 = _mm_sub_ps(C6,V03);

    V00 = XM_PERMUTE_PS(MT.r[1],_MM_SHUFFLE(0,3,0,3));
    // V10 = D0Z,D0Z,D2X,D2Y
    V10 = _mm_shuffle_ps(D0,D2,_MM_SHUFFLE(1,0,2,2));
    V10 = XM_PERMUTE_PS(V10,_MM_SHUFFLE(0,2,3,0));
    V01 = XM_PERMUTE_PS(MT.r[0],_MM_SHUFFLE(2,0,3,1));
    // V11 = D0X,D0W,D2X,D2Y
    V11 = _mm_shuffle_ps(D0,D2,_MM_SHUFFLE(1,0,3,0));
    V11 = XM_PERMUTE_PS(V11,_MM_SHUFFLE(2,1,0,3));
    V02 = XM_PERMUTE_PS(MT.r[3],_MM_SHUFFLE(0,3,0,3));
    // V12 = D1Z,D1Z,D2Z,D2W
    V12 = _mm_shuffle_ps(D1,D2,_MM_SHUFFLE(3,2,2,2));
    V12 = XM_PERMUTE_PS(V12,_MM_SHUFFLE(0,2,3,0));
    V03 = XM_PERMUTE_PS(MT.r[2],_MM_SHUFFLE(2,0,3,1));
    // V13 = D1X,D1W,D2Z,D2W
    V13 = _mm_shuffle_ps(D1,D2,_MM_SHUFFLE(3,2,3,0));
    V13 = XM_PERMUTE_PS(V13,_MM_SHUFFLE(2,1,0,3));

    V00 = _mm_mul_ps(V00,V10);
    V01 = _mm_mul_ps(V01,V11);
    V02 = _mm_mul_ps(V02,V12);
    V03 = _mm_mul_ps(V03,V13);
    XMVECTOR C1 = _mm_sub_ps(C0,V00);
    C0 = _mm_add_ps(C0,V00);
    XMVECTOR C3 = _mm_add_ps(C2,V01);
    C2 = _mm_sub_ps(C2,V01);
    XMVECTOR C5 = _mm_sub_ps(C4,V02);
    C4 = _mm_add_ps(C4,V02);
    XMVECTOR C7 = _mm_add_ps(C6,V03);
    C6 = _mm_sub_ps(C6,V03);

    C0 = _mm_shuffle_ps(C0,C1,_MM_SHUFFLE(3,1,2,0));
    C2 = _mm_shuffle_ps(C2,C3,_MM_SHUFFLE(3,1,2,0));
    C4 = _mm_shuffle_ps(C4,C5,_MM_SHUFFLE(3,1,2,0));
    C6 = _mm_shuffle_ps(C6,C7,_MM_SHUFFLE(3,1,2,0));
    C0 = XM_PERMUTE_PS(C0,_MM_SHUFFLE(3,1,2,0));
    C2 = XM_PERMUTE_PS(C2,_MM_SHUFFLE(3,1,2,0));
    C4 = XM_PERMUTE_PS(C4,_MM_SHUFFLE(3,1,2,0));
    C6 = XM_PERMUTE_PS(C6,_MM_SHUFFLE(3,1,2,0));
    // Get the determinate
    XMVECTOR vTemp = XMVector4Dot(C0,MT.r[0]);
    if (pDeterminant != nullptr)
        *pDeterminant = vTemp;
    //vTemp = _mm_div_ps(g_XMOne,vTemp);
    vTemp = _mm_div_ps(_mm_set1_ps(1.0f),vTemp);
    XMMATRIX mResult;
    mResult.r[0] = _mm_mul_ps(C0,vTemp);
    mResult.r[1] = _mm_mul_ps(C2,vTemp);
    mResult.r[2] = _mm_mul_ps(C4,vTemp);
    mResult.r[3] = _mm_mul_ps(C6,vTemp);
    return mResult;
}

//------------------------------------------------------------------------------

inline XMMATRIX XMMatrixPerspectiveFovRH
(
    float FovAngleY,
    float AspectHByW,
    float NearZ,
    float FarZ
)
{
    float    SinFov;
    float    CosFov;
    XMScalarSinCos(&SinFov, &CosFov, 0.5f * FovAngleY);
    float fRange = FarZ / (NearZ-FarZ);
    // Note: This is recorded on the stack
    float Height = CosFov / SinFov;
    XMVECTOR rMem = {
        Height / AspectHByW,
        Height,
        fRange,
        fRange * NearZ
    };
    // Copy from memory to SSE register
    XMVECTOR vValues = rMem;
    XMVECTOR vTemp = _mm_setzero_ps();
    // Copy x only
    vTemp = _mm_move_ss(vTemp,vValues);
    // CosFov / SinFov,0,0,0
    XMMATRIX M;
    M.r[0] = vTemp;
    // 0,Height / AspectHByW,0,0
    vTemp = vValues;
    vTemp = _mm_and_ps(vTemp,g_XMMaskY);
    M.r[1] = vTemp;
    // x=fRange,y=-fRange * NearZ,0,-1.0f
    vTemp = _mm_setzero_ps();
    //vValues = _mm_shuffle_ps(vValues,g_XMNegIdentityR3,_MM_SHUFFLE(3,2,3,2));
	vValues = _mm_shuffle_ps(vValues,float4(0,0,0,-1.0f).v,_MM_SHUFFLE(3,2,3,2));
    // 0,0,fRange,-1.0f
    vTemp = _mm_shuffle_ps(vTemp,vValues,_MM_SHUFFLE(3,0,0,0));
    M.r[2] = vTemp;
    // 0,0,fRange * NearZ,0.0f
    vTemp = _mm_shuffle_ps(vTemp,vValues,_MM_SHUFFLE(2,1,0,0));
    M.r[3] = vTemp;
    return M;
}

//------------------------------------------------------------------------------

inline XMMATRIX XMMatrixLookToLH
(
    FXMVECTOR EyePosition,
    FXMVECTOR EyeDirection,
    FXMVECTOR UpDirection
)
{
	XMVECTOR R2 = XMVector3Normalize(EyeDirection);

	XMVECTOR R0 = XMVector3Cross(UpDirection, R2);
	R0 = XMVector3Normalize(R0);

	XMVECTOR R1 = XMVector3Cross(R2, R0);

	XMVECTOR NegEyePosition = XMVectorNegate(EyePosition);

	XMVECTOR D0 = XMVector3Dot(R0, NegEyePosition);
	XMVECTOR D1 = XMVector3Dot(R1, NegEyePosition);
	XMVECTOR D2 = XMVector3Dot(R2, NegEyePosition);

	XMMATRIX M;
	M.r[0] = float4::select(float4(D0),float4(R0),float4::selectctrl(1,1,1,0)).v;//XMVectorSelect(D0, R0, g_XMSelect1110.v);
	M.r[1] = float4::select(float4(D1),float4(R1),float4::selectctrl(1,1,1,0)).v;//XMVectorSelect(D1, R1, g_XMSelect1110.v);
	M.r[2] = float4::select(float4(D2),float4(R2),float4::selectctrl(1,1,1,0)).v;//XMVectorSelect(D2, R2, g_XMSelect1110.v);
	M.r[3] = float4(0,0,0,1.0f).v;//g_XMIdentityR3.v;

	M = XMMatrixTranspose(M);

	return M;
}

inline XMMATRIX XMMatrixLookToRH
(
    FXMVECTOR EyePosition,
    FXMVECTOR EyeDirection,
    FXMVECTOR UpDirection
)
{
    XMVECTOR NegEyeDirection = XMVectorNegate(EyeDirection);
    return XMMatrixLookToLH(EyePosition, NegEyeDirection, UpDirection);
}

#undef XMVECTOR
#undef FXMVECTOR

#endif // SMPORT_INL
