/******************************************************************************
 * libKrylov
 * =========
 * NKPF.hxx
 *
 * This file contains definitions pertaining to Newton-Krylov power flow
 *
 * 21 August 2014
 * ~ ry
 * ***************************************************************************/
#ifndef KRY_NKPF_HXX
#define KRY_NKPF_HXX

#include "kry/Math.hxx"

namespace kry
{

class NKPF;
struct jidx;

struct jidx { size_t j0, j1; };
struct JacobiMap
{
  size_t j0_sz, j1_sz;
  std::vector<jidx> map;
  size_t size();
};

class NKPF
{
  public:

    NKPF(SparseMatrix Y, SparseMatrix YA, JacobiMap jmap, size_t n);
    
    Matrix Q, //Subspace Projector
           H; //Hessenburg Reduction

    Vector ve,  //voltage estimate
           dve, //voltage delta estimate
           ps,  //scheduled power
           pc,  //calculated power
           dp,  //power delta
           dv0, //voltage delta initial guess
           dr0, //initial residual of the Jacobean system
           qdp, //projected power delta
           qdv; //projected voltage delta

    SparseMatrix Y,  //admittance matrix magnitudes
                 YA; //admittance matrix angles
    
    JacobiMap jmap; //map bus indices onto Jacobi indices

    //the power flow equations
    double p(size_t), q(size_t); 

    /* power gradient functions
    -----------------------------------------------------*/
    //real-power gradient
    double jdp(size_t), 
           jdp_va(size_t), 
           jdp_v(size_t),
           dp_dva(size_t, size_t), 
           dp_dv(size_t, size_t);
    
    //reactive-power gradient
    double jdq(size_t), 
           jdq_va(size_t), 
           jdq_v(size_t), 
           dq_dva(size_t, size_t), 
           dq_dv(size_t, size_t);
    /*---------------------------------------------------*/

    //accessors for voltage magnitude and angle
    double v(size_t), va(size_t);

    //accessors for voltage delta magnitude and angle
    double dv(size_t), dva(size_t);

    //conductance and susceptance per bus
    double g(size_t), b(size_t);
     
};

}

#endif