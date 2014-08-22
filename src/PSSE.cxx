/******************************************************************************
 * libKrylov
 * =========
 * PSSE.cxx contains the implementation of the parsing and utilization
 * facilities
 *
 * 20 August 2014
 * ~ ry
 * ***************************************************************************/
#include "kry/PSSE.hxx"

using namespace kry;
using namespace kry::input::psse;
using namespace kry::input;
using std::string;
using std::ifstream;
using std::istreambuf_iterator;
using std::stringstream;
using std::endl;
using std::stoul;
using std::stod;
using std::invalid_argument;
using std::cout;
using std::runtime_error;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
template <>
double kry::input::psse::parseElem<double>(const string & src,
                         const string & err,
                         ulong idx)
{ return std::stod(src); }

template <>
unsigned long kry::input::psse::parseElem<unsigned long>(const string & src,
                                       const string & err,
                                       unsigned long idx)
{ return std::stoul(src); }

template <>
unsigned short kry::input::psse::parseElem<unsigned short>(const string & src,
                                        const string & err,
                                        unsigned long idx)
{ return std::stoul(src); }
#pragma clang diagnostic pop

kry::input::psse::Source::Source(string fn)
{
    ifstream ifs = ifstream(fn);
    text = string(istreambuf_iterator<char>(ifs),
                         istreambuf_iterator<char>());
    
    computeLayout();
    resolveTransfomerIndicies();
    parseElements();
    assignSystemIds();
}

void kry::input::psse::Source::parseBus(ulong idx)
{
    ulong number;
    string name;
    double base_kv, vm, va;
    unsigned short type, area, zone, owner;
    complex y{}, v;
    
    string str = text.substr(idx, nextRow(text, idx) - idx);
    vector<string> toks = kry::input::split(str, ',');
    ulong cidx{idx};
    
    //number
    number = parseElem<ulong>(toks[0], "Invalid bus number", cidx);
    cidx += toks[0].length();
    //name
    name = toks[1];
    cidx += toks[1].length();
    //base_kv
    base_kv = parseElem<double>(toks[2], "Invalid base kv", cidx);
    cidx += toks[2].length();
    //type
    type = parseElem<unsigned short>(toks[3], "Invalid bus type", cidx);
    cidx += toks[3].length();
    //area
    area = parseElem<unsigned short>(toks[4], "Invalid bus area", cidx);
    cidx += toks[4].length();
    //zone
    zone = parseElem<unsigned short>(toks[5], "Invalid bus zone", cidx);
    cidx += toks[5].length();
    //owner
    owner = parseElem<unsigned short>(toks[6], "Invalid owner number", cidx);
    cidx += toks[6].length();
    //vm
    vm = parseElem<double>(toks[7], "Invalid bus voltage magnitude", cidx);
    cidx += toks[7].length();
    //va
    va = parseElem<double>(toks[8], "Invalid bus voltage angle", cidx);
    cidx += toks[8].length();
    
    v = std::polar(vm, (va * M_PI/180.0));
    //auto loc = sourceCoordsFromAbsolute(text, idx);
    Bus b{number, base_kv, type, y, v};
    b.name = name;
    b.area = area;
    b.zone = zone;
    b.owner = owner;
    
    buses.push_back(b);
}

void kry::input::psse::Source::parseHeader()
{
  string str = text.substr(0, nextRow(text, 0));
  vector<string> toks = kry::input::split(str, ',');
  base_mva = parseElem<double>(toks[1], "Invalid MVA base", 0); 
}

void kry::input::psse::Source::parseLoad(ulong idx)
{
    ulong number;
    string id;
    complex power;
    unsigned short status, area, zone, owner;
    double p, q, ip, iq, yp, yq;
    
    string str = text.substr(idx, nextRow(text, idx) - idx);
    vector<string> toks = kry::input::split(str, ',');
    ulong cidx{idx};
    
    //number
    number = parseElem<ulong>(toks[0], "Invalid load number", cidx);
    cidx += toks[0].length();
    //id
    id = toks[1];
    cidx += toks[1].length();
    //status
    status = parseElem<unsigned short>(toks[2], "Invalid load status", cidx);
    cidx += toks[2].length();
    //area
    area = parseElem<unsigned short>(toks[3], "Invalid load area", cidx);
    cidx += toks[3].length();
    //zone
    zone = parseElem<unsigned short>(toks[4], "Invalid load zone", cidx);
    cidx += toks[4].length();
    //load power
    p = parseElem<double>(toks[5], "Invalid load real power", cidx);
    cidx += toks[5].length();
    q = parseElem<double>(toks[6], "Invalid load reactive power", cidx);
    cidx += toks[6].length();
    //load constant current power
    ip = parseElem<double>(toks[7],
                           "invalid load real constant-current power",
                           cidx);
    cidx += toks[7].length();
    iq = parseElem<double>(toks[8],
                           "Invalid load reactive constant-current power",
                           cidx);
    cidx += toks[8].length();
    //load constant impedance power
    yp = parseElem<double>(toks[9],
                           "Invalid load real constant-impedance power",
                           cidx);
    cidx += toks[9].length();
    yq = parseElem<double>(toks[10],
                           "Invalid load reactive constant-impedance power",
                           cidx);
    cidx += toks[10].length();
    //owner
    owner = parseElem<unsigned short>(toks[11], "Invalid load owner", cidx);
    cidx += toks[11].length();
    
    Load l{number, complex{p, q}, status};
    l.id = id;
    l.area = area;
    l.zone = zone;
    l.owner = owner;
    l.ipq = complex{ip, iq};
    l.ypq = complex{yp, yq};
    
    loads.push_back(l);
}

void kry::input::psse::Source::parseGen(ulong idx)
{
    ulong number, ireg;
    string id;
    complex p, zz, zt;
    double pg, qg, qt, qb, vs, mbase, gtap, rmpct, pt, pb, zr, zx, rt, xt;
    unsigned short stat;
    
    string str = text.substr(idx, nextRow(text, idx) - idx);
    vector<string> toks = kry::input::split(str, ',');
    ulong cidx{idx};
    
    //number
    number = parseElem<ulong>(toks[0], "Invalid generator number", cidx);
    cidx += toks[0].length();
    //id
    id = toks[1];
    cidx += toks[1].length();
    //gen power
    pg = parseElem<double>(toks[2], "Invalid generator real power", cidx);
    cidx += toks[2].length();
    qg = parseElem<double>(toks[3], "Invalid generator reactive power", cidx);
    cidx += toks[3].length();
    //gen reactive limits
    qt = parseElem<double>(toks[4], "Invalid generator reactive power max",
                           cidx);
    cidx += toks[4].length();
    qb = parseElem<double>(toks[5], "Invalid generator reactive power min",
                           cidx);
    cidx += toks[5].length();
    //regulated voltage setpoint
    vs = parseElem<double>(toks[6],
                           "Invalid generator regulated voltage setpoint",
                           cidx);
    cidx += toks[6].length();
    //regulated bus id
    ireg = parseElem<ulong>(toks[7], "Invalid generator regulated bus id",
                            cidx);
    cidx += toks[7].length();
    //MVA base
    mbase = parseElem<double>(toks[8], "Invalid generator MVA base", cidx);
    cidx += toks[8].length();
    //gen impedance
    zr = parseElem<double>(toks[9], "Invalid generator resistance", cidx);
    cidx += toks[9].length();
    zx = parseElem<double>(toks[10], "Invalid generator reacatnce", cidx);
    cidx += toks[10].length();
    //step up impedance
    rt = parseElem<double>(toks[11], "Invalid generator step-up resistance",
                           cidx);
    cidx += toks[11].length();
    xt = parseElem<double>(toks[12], "Invalid generator step-up reactance",
                           cidx);
    cidx += toks[12].length();
    //gtap
    gtap = parseElem<double>(toks[12], "Invalid generator gtap", cidx);
    cidx += toks[12].length();
    //status
    stat = parseElem<unsigned short>(toks[13], "Invalid generator status",
                                     cidx);
    cidx += toks[13].length();
    //rmpct
    rmpct = parseElem<double>(toks[14], "Invalid generator rmpct", cidx);
    cidx += toks[14].length();
    //PT/PB
    pt = parseElem<double>(toks[15], "Invalid generator max real power", cidx);
    cidx += toks[15].length();
    pb = parseElem<double>(toks[16], "Invalid generator min real power", cidx);
    cidx += toks[16].length();
    
    Gen g{number, complex{pg, qg}, stat};
    g.ireg = ireg;
    g.id = id;
    g.qt = qt;
    g.qb = qb;
    g.vs = vs;
    g.mbase = mbase;
    g.gtap = gtap;
    g.rmpct = rmpct;
    g.pt = pt;
    g.pb = pb;
    
    gens.push_back(g);
}

void kry::input::psse::Source::parseLine(ulong idx)
{
    ulong i, j;
    string ckt;
    double b, ratea, rateb, ratec, len, zr, zx, yig, yib, yjg, yjb;
    unsigned short st;
    
    string str = text.substr(idx, nextRow(text, idx) - idx);
    vector<string> toks = kry::input::split(str, ',');
    ulong cidx{idx};
    
    //bus id's
    i = parseElem<ulong>(toks[0], "Invalid line ibus id", cidx);
    cidx += toks[0].length();
    j = parseElem<ulong>(toks[1], "Invalid line jbus id", cidx);
    cidx += toks[1].length();
    //circuit identifier
    ckt = toks[2];
    cidx += toks[2].length();
    //line impedance
    zr = parseElem<double>(toks[3], "Invalid line resistance", cidx);
    cidx += toks[3].length();
    zx = parseElem<double>(toks[4], "Invalid line reactance", cidx);
    cidx += toks[4].length();
    //line charging susceptance
    b = parseElem<double>(toks[5], "Invalid line charging susceptance", cidx);
    cidx += toks[5].length();
    //line phase ratings
    ratea = parseElem<double>(toks[6], "Invalid a-phase line rating", cidx);
    cidx += toks[6].length();
    rateb = parseElem<double>(toks[7], "Invalid b-phase line rating", cidx);
    cidx += toks[7].length();
    ratec = parseElem<double>(toks[8], "Invalid c-phase line rating", cidx);
    cidx += toks[8].length();
    //line terminal shunt admittances
    //i
    yig = parseElem<double>(toks[9],
                            "Invalid line i-terminal shunt conductance", cidx);
    cidx += toks[9].length();
    yib = parseElem<double>(toks[10],
                            "Invalid line i-terminal shunt susceptance", cidx);
    cidx += toks[10].length();
    //j
    yjg = parseElem<double>(toks[11],
                            "Invalid line j-terminal shunt conductance", cidx);
    cidx += toks[11].length();
    yjb = parseElem<double>(toks[12],
                            "Invalid line j-terminal shunt susceptance", cidx);
    cidx += toks[12].length();
    //line status
    st = parseElem<unsigned short>(toks[13], "Invalid line status", cidx);
    cidx += toks[13].length();
    //line length
    len = parseElem<double>(toks[14], "Invalid line length", cidx);
    cidx += toks[14].length();
    
    Line l{i, j, complex{zr, zx}};
    l.ckt = ckt;
    l.b = b;
    l.ratea = ratea;
    l.rateb = rateb;
    l.ratec = ratec;
    l.len = len;
    l.yi = complex{yig, yib};
    l.yj = complex{yjg, yjb};
    l.st = st;
    
    lines.push_back(l);
}

void kry::input::psse::Source::parseShunt(ulong idx)
{
  FixedShunt fs;

  string str = text.substr(idx, nextRow(text, idx) - idx);
  vector<string> toks = kry::input::split(str, ',');
  ulong cidx{idx};

  fs.i = parseElem<ulong>(toks[0], "Invalid Fixed shunt bus", cidx);
  cidx += toks[0].length();

  cidx += toks[1].length();
  cidx += toks[2].length();
  cidx += toks[3].length();

  fs.b = parseElem<ulong>(toks[4], "Invalid Fixed shunt susceptance", cidx);
  fs.b /= base_mva;

  shunts.push_back(fs);
}

void kry::input::psse::Source::parseTransformer(ulong idx)
{
    unsigned short type;
    
    string str = text.substr(idx, nextRow(text, idx) - idx);
    vector<string> toks = kry::input::split(str, ',');
    ulong cidx{idx};
    
    //record1
    ulong i, j, k;
    string ckt, name;
    unsigned short cw, cz, cm, nmetr, stat;
    double magnetizing_g, magnetizing_b;
    
    //i,j,k bus
    i = parseElem<ulong>(toks[0], "Invalid transformer i-bus", cidx);
    cidx += toks[0].length();
    j = parseElem<ulong>(toks[1], "Invalid transformer j-bus", cidx);
    cidx += toks[1].length();
    k = parseElem<ulong>(toks[2], "Invalid transformer k-bus", cidx);
    cidx += toks[2].length();
    
    //ckt id
    ckt = toks[3];
    cidx += toks[3].length();
    
    //cw, cz, cm
    cw = parseElem<unsigned short>(toks[4], "Invalid transformer cw code", cidx);
    cidx += toks[4].length();
    cz = parseElem<unsigned short>(toks[5], "Invalid transformer cz code", cidx);
    cidx += toks[5].length();
    cm = parseElem<unsigned short>(toks[6], "Invalid transformer cm code", cidx);
    cidx += toks[6].length();
    
    //magnetizing conductance/susceptance
    magnetizing_g = parseElem<double>(toks[7],
                                "Invalid transformer magnetizing conductance",
                                      cidx);
    cidx += toks[7].length();
    magnetizing_b = parseElem<double>(toks[8],
                                "Invalid transformer magnetizing susceptance",
                                      cidx);
    cidx += toks[8].length();
    
    //nmetr
    nmetr = parseElem<unsigned short>(toks[9],
                                      "Invalid transformer nmeter", cidx);
    cidx += toks[9].length();
    
    //name
    name = toks[10];
    
    stat = parseElem<unsigned short>(toks[11], "Invalid transformer status",
                                     cidx);
    cidx += toks[11].length();
    
    
    //record2
    double z12r, z12x, z23r, z23x, z31r, z31x;
    double sbase12, sbase23, sbase31, vmstar, anstar;
    
    idx = nextRow(text, idx);
    str = text.substr(idx, nextRow(text, idx) - idx);
    toks.clear();
    toks = kry::input::split(str, ',');
    cidx = idx;
    
    //z12
    z12r = parseElem<double>(toks[0],
                             "Invalid transformer 12 winding resistance", cidx);
    cidx += toks[0].length();
    z12x = parseElem<double>(toks[1],
                             "Invalid transformer 12 winding reactance", cidx);
    cidx += toks[1].length();
    
    //sbase12
    sbase12 = parseElem<double>(toks[2],
                             "Invalid transformer 12 base MVA", cidx);
    cidx += toks[2].length();
    
    type = toks.size() > 3 ? 5 : 4;
    
    if(type == 5)
    {
        //z23
        z23r = parseElem<double>(toks[3],
                                 "Invalid transformer 23 winding resistance",
                                 cidx);
        cidx += toks[3].length();
        z23x = parseElem<double>(toks[4],
                                 "Invalid transformer 23 winding reactance",
                                 cidx);
        cidx += toks[4].length();
        
        //sbase23
        sbase23 = parseElem<double>(toks[5],
                                    "Invalid transformer 23 base MVA", cidx);
        cidx += toks[5].length();
        
        //z31
        z31r = parseElem<double>(toks[6],
                                 "Invalid transformer 31 winding resistance",
                                 cidx);
        cidx += toks[6].length();
        z31x = parseElem<double>(toks[7],
                                 "Invalid transformer 31 winding reactance",
                                 cidx);
        cidx += toks[7].length();
        
        //sbase23
        sbase31 = parseElem<double>(toks[8],
                                    "Invalid transformer 31 base MVA", cidx);
        cidx += toks[8].length();
        
        vmstar = parseElem<double>(toks[9], "Invalid transformer vmstar", cidx);
        cidx += toks[9].length();
        
        anstar = parseElem<double>(toks[10], "Invalid transformer anstar", cidx);
        cidx += toks[10].length();
    }
    
    //record3
    double windv1, nomv1, ang1, rata1, ratb1, ratc1, rma1, rmi1, vma1, vmi1;
    unsigned short cod1, ntp1, tab1;
    ulong cont1;
    double cx1, cr1;
    
    idx = nextRow(text, idx);
    str = text.substr(idx, nextRow(text, idx) - idx);
    toks.clear();
    toks = kry::input::split(str, ',');
    cidx = idx;
    
    windv1 = parseElem<double>(toks[0], "Invalid transformer windv1", cidx);
    cidx += toks[0].length();
    
    nomv1 = parseElem<double>(toks[1], "Invalid transformer nomv1", cidx);
    cidx += toks[1].length();
    
    ang1 = parseElem<double>(toks[2], "Invalid transformer ang1", cidx);
    cidx += toks[2].length();
    
    rata1 = parseElem<double>(toks[3], "Invalid transformer rata1", cidx);
    cidx += toks[3].length();
    
    ratb1 = parseElem<double>(toks[4], "Invalid transformer ratb1", cidx);
    cidx += toks[4].length();
    
    ratc1 = parseElem<double>(toks[5], "Invalid transformer ratc1", cidx);
    cidx += toks[5].length();
    
    cod1 = parseElem<unsigned short>(toks[6], "Invalid transformer cod1", cidx);
    cidx += toks[6].length();
    
    cont1 = parseElem<ulong>(toks[7], "Invalid transformer cont1", cidx);
    cidx += toks[7].length();
    
    rma1 = parseElem<double>(toks[8], "Invalid transformer rma1", cidx);
    cidx += toks[8].length();
    
    rmi1 = parseElem<double>(toks[9], "Invalid transformer rmi1", cidx);
    cidx += toks[9].length();
    
    vma1 = parseElem<double>(toks[10], "Invalid transformer vma1", cidx);
    cidx += toks[10].length();
    
    vmi1 = parseElem<double>(toks[11], "Invalid transformer vmi1", cidx);
    cidx += toks[11].length();
    
    ntp1 = parseElem<unsigned short>(toks[12], "Invalid transformer ntp1", cidx);
    cidx += toks[12].length();
    
    tab1 = parseElem<unsigned short>(toks[13], "Invalid transformer tab1", cidx);
    cidx += toks[13].length();
    
    cr1 = parseElem<double>(toks[14], "Invalid transformer cr1", cidx);
    cidx += toks[14].length();
    
    cx1 = parseElem<double>(toks[15], "Invalid transformer cx1", cidx);
    cidx += toks[14].length();
    
    
    
    //record4
    double windv2, nomv2, ang2, rata2, ratb2, ratc2, rma2, rmi2, vma2, vmi2;
    unsigned short cod2, ntp2, tab2;
    ulong cont2;
    double cr2, cx2;
    
    idx = nextRow(text, idx);
    str = text.substr(idx, nextRow(text, idx) - idx);
    toks.clear();
    toks = kry::input::split(str, ',');
    cidx = idx;
    
    windv2 = parseElem<double>(toks[0], "Invalid transformer windv2", cidx);
    cidx += toks[0].length();
    
    nomv2 = parseElem<double>(toks[1], "Invalid transformer nomv2", cidx);
    cidx += toks[1].length();
    
    if(type == 5)
    {
    
        ang2 = parseElem<double>(toks[2], "Invalid transformer ang2", cidx);
        cidx += toks[2].length();
        
        rata2 = parseElem<double>(toks[3], "Invalid transformer rata2", cidx);
        cidx += toks[3].length();
        
        ratb2 = parseElem<double>(toks[4], "Invalid transformer ratb2", cidx);
        cidx += toks[4].length();
        
        ratc2 = parseElem<double>(toks[5], "Invalid transformer ratc2", cidx);
        cidx += toks[5].length();
        
        cod2 = parseElem<unsigned short>(toks[6],
                                         "Invalid transformer cod2", cidx);
        cidx += toks[6].length();
        
        cont2 = parseElem<ulong>(toks[7], "Invalid transformer cont2", cidx);
        cidx += toks[7].length();
        
        rma2 = parseElem<double>(toks[8], "Invalid transformer rma2", cidx);
        cidx += toks[8].length();
        
        rmi2 = parseElem<double>(toks[9], "Invalid transformer rmi2", cidx);
        cidx += toks[9].length();
        
        vma2 = parseElem<double>(toks[10], "Invalid transformer vma2", cidx);
        cidx += toks[10].length();
        
        vmi2 = parseElem<double>(toks[11], "Invalid transformer vmi2", cidx);
        cidx += toks[11].length();
        
        ntp2 = parseElem<unsigned short>(toks[12],
                                         "Invalid transformer ntp2", cidx);
        cidx += toks[12].length();
        
        tab2 = parseElem<unsigned short>(toks[13],
                                         "Invalid transformer tab2", cidx);
        cidx += toks[13].length();
        
        cr2 = parseElem<double>(toks[14], "Invalid transformer cr2", cidx);
        cidx += toks[14].length();
        
        cx2 = parseElem<double>(toks[15], "Invalid transformer cx2", cidx);
        cidx += toks[14].length();
    }

    //record5
    double windv3, nomv3, ang3, rata3, ratb3, ratc3, rma3, rmi3, vma3, vmi3;
    unsigned short cod3, ntp3, tab3;
    ulong cont3;
    double cr3, cx3;
    
    if(type == 5)
    {
        
        
        idx = nextRow(text, idx);
        str = text.substr(idx, nextRow(text, idx) - idx);
        toks.clear();
        toks = kry::input::split(str, ',');
        cidx = idx;
        
        windv3 = parseElem<double>(toks[0], "Invalid transformer windv3", cidx);
        cidx += toks[0].length();
        
        nomv3 = parseElem<double>(toks[1], "Invalid transformer nomv3", cidx);
        cidx += toks[1].length();
        
        ang3 = parseElem<double>(toks[2], "Invalid transformer ang3", cidx);
        cidx += toks[2].length();
        
        rata3 = parseElem<double>(toks[3], "Invalid transformer rata3", cidx);
        cidx += toks[3].length();
        
        ratb3 = parseElem<double>(toks[4], "Invalid transformer ratb3", cidx);
        cidx += toks[4].length();
        
        ratc3 = parseElem<double>(toks[5], "Invalid transformer ratc3", cidx);
        cidx += toks[5].length();
        
        cod3 = parseElem<unsigned short>(toks[6],
                                         "Invalid transformer cod3", cidx);
        cidx += toks[6].length();
        
        cont3 = parseElem<ulong>(toks[7], "Invalid transformer cont3", cidx);
        cidx += toks[7].length();
        
        rma3 = parseElem<double>(toks[8], "Invalid transformer rma3", cidx);
        cidx += toks[8].length();
        
        rmi3 = parseElem<double>(toks[9], "Invalid transformer rmi3", cidx);
        cidx += toks[9].length();
        
        vma3 = parseElem<double>(toks[10], "Invalid transformer vma3", cidx);
        cidx += toks[10].length();
        
        vmi3 = parseElem<double>(toks[11], "Invalid transformer vmi3", cidx);
        cidx += toks[11].length();
        
        ntp3 = parseElem<unsigned short>(toks[12], "Invalid transformer ntp3", cidx);
        cidx += toks[12].length();
        
        tab3 = parseElem<unsigned short>(toks[13], "Invalid transformer tab3", cidx);
        cidx += toks[13].length();
        
        cr3 = parseElem<double>(toks[14], "Invalid transformer cr3", cidx);
        cidx += toks[14].length();
        
        cx3 = parseElem<double>(toks[15], "Invalid transformer cx3", cidx);
        cidx += toks[14].length();
    }
    
    Transformer t{i, j, complex{z12r, z12x}, windv1, stat};
    t.type = type;
    //record 1
    t.k = k;
    t.ckt = ckt;
    t.cw = cw;
    t.cz = cz;
    t.cm = cm;
    t.magnetizing_y = complex{magnetizing_g, magnetizing_b};
    t.nmetr = nmetr;
    t.name = name;
    
    //record 2
    t.z12 = complex{z12r, z12x};
    t.sbase12 = sbase12;
    if(type == 5)
    {
        t.z23 = complex{z23r, z23x};
        t.sbase23 = sbase23;
        t.z31 = complex{z31r, z31x};
        t.sbase31 = sbase31;
        t.vmstar = vmstar;
        t.anstar = anstar;
    }
    
    //record 3
    t.windv1 = windv1;
    t.nomv1 = nomv1;
    t.ang1 = ang1;
    t.rata1 = rata1;
    t.ratb1 = ratb1;
    t.ratc1 = ratc1;
    t.cod1 = cod1;
    t.cont1 = cont1;
    t.rma1  = rma1;
    t.rmi1 = rmi1;
    t.vma1 = vma1;
    t.vmi1 = vmi1;
    t.ntp1 = ntp1;
    t.tab1 = tab1;
    t.cz1 = complex{cr1, cx1};
    
    //record 4
    t.windv2 = windv2;
    t.nomv2 = nomv2;
    if(type == 5)
    {
        t.ang2 = ang2;
        t.rata2 = rata2;
        t.ratb2 = ratb2;
        t.ratc2 = ratc2;
        t.cod2 = cod2;
        t.cont2 = cont2;
        t.rma2  = rma2;
        t.rmi2 = rmi2;
        t.vma2 = vma2;
        t.vmi2 = vmi2;
        t.ntp2 = ntp2;
        t.tab2 = tab2;
        t.cz2 = complex{cr2, cx2};
    }
    
    //record 5
    if(type == 5)
    {
        t.windv3 = windv3;
        t.nomv3 = nomv3;
        t.ang3 = ang3;
        t.rata3 = rata3;
        t.ratb3 = ratb3;
        t.ratc3 = ratc3;
        t.cod3 = cod3;
        t.cont3 = cont3;
        t.rma3  = rma3;
        t.rmi3 = rmi3;
        t.vma3 = vma3;
        t.vmi3 = vmi3;
        t.ntp3 = ntp3;
        t.tab3 = tab3;
        t.cz3 = complex{cr3, cx3};
    }
    
    transformers.push_back(t);
}

void kry::input::psse::Source::parseElements()
{
    for(auto r : bus_rows){ parseBus(r); }
    std::sort(buses.begin(), buses.end(),
              [](const Bus &a, const Bus &b) -> bool{
                  return a.number < b.number;
              }
             );
   
    parseHeader();
    for(auto r : load_rows) { parseLoad(r); }
    for(auto r : gen_rows) { parseGen(r); }
    for(auto r : shunt_rows) { parseShunt(r); }
    for(auto r : line_rows) { parseLine(r); }
    for(auto r : tfmr_rows) { parseTransformer(r); }
}

void kry::input::psse::Source::computeLayout()
{
    ulong
    buses_begin = rowIndex(text, 3),
    buses_end = endOfSection(text, buses_begin, "End of Bus data"),
    loads_begin = startOfSection(text, buses_end, "Begin Load data"),
    loads_end = endOfSection(text, loads_begin, "End of Load data"),
    shunt_begin = startOfSection(text, loads_end, "Begin Fixed Shunt data"),
    shunt_end = endOfSection(text, shunt_begin, "End of Fixed Shunt data"),
    gens_begin = startOfSection(text, shunt_end, "Begin Generator data"),
    gens_end = endOfSection(text, gens_begin, "End of Generator data"),
    lines_begin = startOfSection(text, gens_end, "Begin Branch data"),
    lines_end = endOfSection(text, lines_begin, "End of Branch data"),
    transformers_begin = startOfSection(text, lines_end,
                                        "Begin Transformer data"),
    transformers_end = endOfSection(text, transformers_begin,
                                    "End of Transformer data");
    
    bus_rows = rowIndices(text, buses_begin, buses_end);
    load_rows = rowIndices(text, loads_begin, loads_end);
    gen_rows = rowIndices(text, gens_begin, gens_end);
    line_rows = rowIndices(text, lines_begin, lines_end);
    tfmr_rows = rowIndices(text, transformers_begin, transformers_end);
    shunt_rows = rowIndices(text, shunt_begin, shunt_end);
}

void kry::input::psse::Source::resolveTransfomerIndicies()
{
    vector<ulong> tfmr{};
    for(ulong i=0; i<tfmr_rows.size(); )
    {
        tfmr.push_back(tfmr_rows[i]);
        if(csvRowSize(text, tfmr_rows[i+1]) > 3) { i += 5; }
        else { i += 4; }
    }
    tfmr_rows = tfmr;
}

string kry::input::psse::Source::basicReport()
{
    stringstream ss{};
    ss << "PsseSource Basic Report ---" << endl
    << "Buses           " << bus_rows.size() << endl
    << "Loads           " << load_rows.size() << endl
    << "Shunts          " << shunt_rows.size() << endl
    << "Generators      " << gen_rows.size() << endl
    << "Lines           " << line_rows.size() << endl
    << "Transformers    " << tfmr_rows.size() << endl;
    return ss.str();
}

Bus & kry::input::psse::Source::findBus(ulong number)
{
    Bus &b = *(std::lower_bound(buses.begin(), buses.end(), number,
                                [](const Bus &b, ulong k) -> bool {
                                    return b.number < k;
                                }
              ));
    
    if(b.number != number)
    {
        stringstream ss{};
        ss << "Phantom bus requested " << number;
        throw runtime_error{ss.str()};
    }
    
    return b;
}

void kry::input::psse::Source::assignSystemIds()
{
  uint seq{0};
  for(auto &b : buses)
  {
    b.sid = seq++;
  }

  for(Line &l : lines)
  {
    l.si = findBus(l.i).sid;
    l.sj = findBus(l.j).sid;
  }
  
  for(Transformer &t : transformers)
  {
    t.si = findBus(t.i).sid;
    t.sj = findBus(t.j).sid;
  }

  for(FixedShunt &s : shunts)
  {
    s.si = findBus(s.i).sid;
  }
}

/*
GridSpace::Grid kry::input::psse::Source::toGrid()
{
  vector<GridSpace::Bus> g_buses{};
  g_buses.reserve(buses.size());
  vector<GridSpace::Load> g_loads{};
  g_loads.reserve(loads.size());
  vector<GridSpace::Generator> g_generators{};
  g_generators.reserve(gens.size());
  vector<GridSpace::Line> g_lines{};
  g_lines.reserve(lines.size());
  vector<GridSpace::Transformer> g_transformers{};
  g_transformers.reserve(transformers.size());
  
  ulong seq{0};
  for(auto &b : buses)
  {
    GridSpace::BusType bt;
    switch(b.type)
    {
        case 1 : bt = GridSpace::BusType::Load; break;
        case 2 : bt = GridSpace::BusType::Generator; break;
        case 3 : bt = GridSpace::BusType::Slack; break;
        case 4 : bt = GridSpace::BusType::Isolated; break;
        default:
        {
            stringstream ss{};
            ss << "Unknown bus type `" << b.type << "` for bus with number "
            << b.number;
            throw runtime_error{ss.str()};
        }
    }
    
    GridSpace::Bus bus{bt, seq};
    b.sid = seq++;
    bus.voltage = b.v;
    bus.shunt_y = b.y;
    bus.name = b.name;
    bus.bus_type = bt;
    g_buses.push_back(bus);
  }
  
  seq = 0;
  for(auto &l : loads)
  {
      Bus &b = findBus(l.number);
      GridSpace::Load load{seq++, b.sid};
      load.power = l.pql;
      g_loads.push_back(load);
  }
  
  seq = 0;
  for(auto &gn : gens)
  {
      Bus &b = findBus(gn.number);
      GridSpace::Generator gen{seq++, b.sid};
      gen.power = gn.p;
      g_generators.push_back(gen);
  }
  
  seq = 0;
  for(auto &l : lines)
  {
      Bus &i = findBus(l.i),
          &j = findBus(l.j);
      GridSpace::Line line{i.sid, j.sid, l.z, l.b};
      line.shunt_i = l.yi;
      line.shunt_j = l.yj;
      g_lines.push_back(line);
      
  }
  
  seq = 0;
  for(auto &t : transformers)
  {
      Bus &i = findBus(t.i),
          &j = findBus(t.j);
      
      double tr = t.windv2 / t.windv1;
      GridSpace::Transformer transformer{i.sid, j.sid, t.z12, tr};
      g_transformers.push_back(transformer);

  }
  
  GridSpace::Grid
  g{std::move(g_buses),
    std::move(g_loads),
    std::move(g_generators),
    std::move(g_lines),
    std::move(g_transformers)};
  
  return std::move(g);
}
*/

std::vector<double> kry::input::psse::Source::staticBusVoltageMagnitudes()
{
  std::vector<double> result; 

  std::transform(buses.begin(), buses.end(), std::back_inserter(result),
      [](const Bus &b)
      {
        return std::abs(b.v);
      });

  return result;
}

std::vector<double> kry::input::psse::Source::staticBusVoltageAngles()
{
  std::vector<double> result; 

  std::transform(buses.begin(), buses.end(), std::back_inserter(result),
      [](const Bus &b)
      {
        return std::arg(b.v);
      });

  return result;
}

std::vector<double> kry::input::psse::Source::pSch()
{
  std::vector<double> result(buses.size(), 0);

  for(auto &l : loads) 
  { 
    l.sid = findBus(l.number).sid;
    result[l.sid] -= l.pql.real(); 
  }

  for(auto &g : gens) 
  { 
    g.sid = findBus(g.number).sid;
    result[g.sid] += g.p.real(); 
  }

  return result;
}

std::vector<double> kry::input::psse::Source::qSch()
{
  std::vector<double> result(buses.size(), 0);

  for(auto &l : loads) 
  { 
    l.sid = findBus(l.number).sid;
    result[l.sid] -= l.pql.imag(); 
  }

  for(auto &g : gens) 
  { 
    g.sid = findBus(g.number).sid;
    result[g.sid] += g.p.imag(); 
  }

  return result;
}

size_t compute_max_node_valence(const Source &s)
{
  vector<size_t> valences(s.buses.size(), 1);
  for(const Line &l : s.lines)
  {
    valences[l.si]++;
    valences[l.sj]++;
  }
  
  for(const Transformer &t : s.transformers)
  {
    valences[t.si]++;
    valences[t.sj]++;
  }

  return *std::max_element(valences.begin(), valences.end());

}

std::array<SparseMatrix, 2> kry::input::psse::ymatrix(const Source &s)
{
  size_t z = compute_max_node_valence(s); 
  size_t n = s.buses.size();

  SparseMatrix  Y(n, n, z),
               YA(n, n, z);

  for(const Line &l : s.lines)
  {
    complex y_i_diag = std::polar(Y(l.si,l.si), YA(l.si,l.si));
    complex y_j_diag = std::polar(Y(l.sj,l.sj), YA(l.sj,l.sj));
    complex y = 1.0/l.z;
    y_i_diag += y + complex(0, 0.5*l.b);
    y_j_diag += y + complex(0, 0.5*l.b);
    Y(l.si, l.si) = abs(y_i_diag);
    Y(l.sj, l.sj) = abs(y_j_diag);
    Y(l.si, l.sj) = Y(l.sj, l.si) = -abs(y);

    YA(l.si, l.si) = arg(y_i_diag);
    YA(l.sj, l.sj) = arg(y_j_diag);
    YA(l.si, l.sj) = YA(l.sj, l.si) = -arg(y);
  }

  for(const Transformer &t : s.transformers)
  {
    complex y_i_diag = std::polar(Y(t.si,t.si), YA(t.si,t.si));
    complex y_j_diag = std::polar(Y(t.sj,t.sj), YA(t.sj,t.sj));
    complex y = 1.0/t.z12;
    double tr = t.windv2 / t.windv1;
    y_i_diag += pow(tr,2) * y;
    y_j_diag += y;

    Y(t.si, t.si) = abs(y_i_diag);
    Y(t.sj, t.sj) = abs(y_j_diag);
    Y(t.si, t.sj) = -tr * abs(y);
    Y(t.sj, t.si) = -tr * abs(y);
    
    YA(t.si, t.si) = arg(y_i_diag);
    YA(t.sj, t.sj) = arg(y_j_diag);
    YA(t.si, t.sj) = -tr * arg(y);
    YA(t.sj, t.si) = -tr * arg(y);
  }

  for(const FixedShunt &s : s.shunts)
  {
    complex y = std::polar(Y(s.si,s.si), YA(s.si,s.si));
    std::cout << "sb" << s.b << std::endl;
    y += complex(0, s.b);
    Y(s.si,s.si) = abs(y);
    YA(s.si,s.si) = arg(y);
  }

  return {{Y, YA}};
}
