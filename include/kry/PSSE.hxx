/******************************************************************************
 * libKrylov
 * =========
 * PSSE.hxx contains data structures for parsing and using PSSE data
 *
 * 20 August 2014
 * ~ ry
 * ***************************************************************************/
#ifndef KRY_ARNOLDI_HXX
#define KRY_ARNOLDI_HXX

#include "Common.hh"
#include "kry/Math.hxx"
#include "NKPF.hxx"

#include <sstream>
#include <iostream>
#include <fstream>
#include <complex>
#include <memory>
#include <array>

namespace kry { namespace input { namespace psse {

using complex = std::complex<double>;
using std::string;
using std::vector;

struct Bus
{
    Bus(ulong number,
        double base_kv,
        unsigned short type,
        complex y,
        complex v)
    : number{number}, base_kv{base_kv}, type{type}, y{y}, v{v} {}
    
    ulong number, sid{}, line_no;
    size_t dva_idx, dv_idx;
    string name{};
    double base_kv;
    unsigned short type, area{}, zone{}, owner{};
    complex y, v;
};

struct Load
{
    Load(ulong number, complex power, unsigned short status)
    : number{number}, status{status}, pql{power} {}
    
    ulong number, sid{};
    string id{};
    unsigned short status, area{}, zone{}, owner{};
    complex pql, ipq{}, ypq{};
};

struct Gen
{
    Gen(ulong number, complex power, unsigned short status)
    : number{number}, status{status}, p{power} {}
    
    ulong number, ireg{}, sid{};
    string id{};
    double qt{}, qb{}, vs{}, mbase{}, gtap{}, rmpct{}, pt{}, pb{}, fi{};
    unsigned short status, oi{};
    complex p, zz{}, zt{};
};

struct Line
{
    Line(ulong i, ulong j, complex z)
    : i{i}, j{j}, z{z} {}
    
    ulong i, j, sid{};
    ulong si, sj;
    string ckt{};
    complex z, yi{}, yj{};
    double b{}, ratea{}, rateb{}, ratec{}, len{}, fi{};
    unsigned short st{}, oi{};
};

struct Transformer
{
    Transformer(ulong i, ulong j, complex z12, double windv1,
                unsigned short status)
    : i{i}, j{j}, stat{status}, z12{z12}, windv1{windv1} {}
    
    unsigned short type;
    
    //record1
    ulong i, j, k{}, sid{};
    ulong si, sj;
    string ckt{}, name{};
    unsigned short cw{}, cz{}, cm{}, nmetr{}, stat{};
    complex magnetizing_y{};
    
    //record2
    complex z12, z23{}, z31{};
    double sbase12{}, sbase23{}, sbase31{}, vmstar{}, anstar{};
    
    //record3
    double windv1, nomv1{}, ang1{}, rata1{}, ratb1{}, ratc1{}, rma1{},
    rmi1{}, vma1{}, vmi1{};
    unsigned short cod1{}, ntp1{}, tab1{};
    ulong cont1{};
    complex cz1{};
    
    //record4
    double windv2{}, nomv2{}, ang2{}, rata2{}, ratb2{}, ratc2{}, rma2{},
    rmi2{}, vma2{}, vmi2{};
    unsigned short cod2{}, ntp2{}, tab2{};
    ulong cont2{};
    complex cz2{};
    
    //record5
    double windv3{}, nomv3{}, ang3{}, rata3{}, ratb3{}, ratc3{}, rma3{},
    rmi3{}, vma3{}, vmi3{};
    unsigned short cod3{}, ntp3{}, tab3{};
    ulong cont3{};
    complex cz3{};
};

struct FixedShunt
{
  size_t i, si;
  double b;
};


template <typename T>
T parseElem(const string & src,
            const string & err,
            unsigned long idx);


struct Source
{
    Source(string fn);

    double base_mva;
    
    vector<ulong>
    bus_rows{}, load_rows{}, gen_rows{}, line_rows{}, tfmr_rows{},
    shunt_rows{};
    
    vector<Bus> buses{};
    vector<Load> loads{};
    vector<Gen> gens{};
    vector<Line> lines{};
    vector<Transformer> transformers{};
    vector<FixedShunt> shunts{};
    
    string text;
    
    string basicReport();
    Bus & findBus(ulong number);

    void assignSystemIds();
    
    //Grid toGrid();
    std::vector<double> staticBusVoltageMagnitudes();
    std::vector<double> staticBusVoltageAngles();
    std::vector<double> pSch();
    std::vector<double> qSch();

    JacobiMap jmap();

    //std::unique_ptr<Model::Protocol::PGrid> getPGrid();
    //std::unique_ptr<Model::Protocol::DGrid> getDGrid();
    
private:
    void computeLayout();
    void resolveTransfomerIndicies();
    void parseElements();
    void parseHeader();
    void parseBus(ulong idx);
    void parseLoad(ulong idx);
    void parseGen(ulong idx);
    void parseLine(ulong idx);
    void parseShunt(ulong idx);
    void parseTransformer(ulong idx);
    
    
    template <typename T>
    T parseElem(const string & src,
                const string & err,
                ulong idx)
    {
        try{ return ::kry::input::psse::parseElem<T>(src, err, idx); }
        catch(std::invalid_argument){
            std::stringstream ss{};
            ss << sourceCoordsFromAbsolute(text, idx)
            << err << " `" << src << "`";
            throw std::runtime_error{ss.str()};
        }
    }
    
};

std::array<SparseMatrix,2> ymatrix(const Source &);


}}} //namespace kry::input::psse

#endif
