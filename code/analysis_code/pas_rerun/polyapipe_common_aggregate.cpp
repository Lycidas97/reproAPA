// Stream accepted cell-by-PAS counts onto the frozen metacell axis.
#include <zlib.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using Index = std::unordered_map<std::string, uint32_t>;
Index index_file(const char* path) {
    Index result; std::ifstream f(path); std::string line;
    if (!f) throw std::runtime_error("Cannot read axis");
    while (std::getline(f,line)) {
        auto tab=line.find('\t');
        if(tab==std::string::npos || !result.emplace(line.substr(0,tab),std::stoul(line.substr(tab+1))).second)
            throw std::runtime_error("Duplicate or malformed axis");
    }
    return result;
}
int main(int argc,char** argv) {
    try {
        if(argc!=8) throw std::runtime_error("counts PAS-map CB-map rows columns output audit");
        auto pas=index_file(argv[2]), cells=index_file(argv[3]);
        uint64_t nr=std::stoull(argv[4]), nc=std::stoull(argv[5]);
        if(nr>UINT32_MAX || nc>UINT32_MAX || nr>SIZE_MAX/nc) throw std::runtime_error("Invalid shape");
        std::vector<uint64_t> matrix(nr*nc,0); std::unordered_map<std::string,bool> observed;
        gzFile f=gzopen(argv[1],"rb"); if(!f) throw std::runtime_error("Cannot open counts");
        gzbuffer(f,1<<20); char line[16384];
        if(!gzgets(f,line,sizeof(line)) || std::string(line)!="gene\tcell\tcount\n") throw std::runtime_error("Count header mismatch");
        uint64_t records=0,total=0;
        while(gzgets(f,line,sizeof(line))) {
            std::string s(line); if(s.back()!='\n') throw std::runtime_error("Truncated or oversized record");
            s.pop_back(); auto a=s.find('\t'),b=s.find('\t',a+1);
            if(a==std::string::npos || b==std::string::npos || s.find('\t',b+1)!=std::string::npos) throw std::runtime_error("Malformed record");
            auto p=pas.find(s.substr(0,a)), c=cells.find(s.substr(a+1,b-a-1));
            if(p==pas.end() || c==cells.end()) throw std::runtime_error("Unknown PAS or frozen barcode");
            auto number=s.substr(b+1); if(number.empty() || number.find_first_not_of("0123456789")!=std::string::npos) throw std::runtime_error("Non-integer count");
            uint64_t value=std::stoull(number),limit=INT64_MAX;
            if(!value || value>limit || c->second>=nr || p->second>=nc) throw std::runtime_error("Invalid count or index");
            auto &slot=matrix[c->second*nc+p->second];
            if(slot>limit-value || total>limit-value) throw std::runtime_error("Integer overflow");
            slot+=value;total+=value;records++;observed[c->first]=true;
        }
        int err=0; gzerror(f,&err); if(err!=Z_OK && err!=Z_STREAM_END) throw std::runtime_error("Invalid gzip");
        if(gzclose(f)!=Z_OK) throw std::runtime_error("Incomplete gzip");
        std::ofstream out(argv[6],std::ios::binary); if(!out) throw std::runtime_error("Cannot write matrix");
        uint64_t nnz=0;
        for(uint32_t r=0;r<nr;r++) for(uint32_t c=0;c<nc;c++) {
            auto value=matrix[uint64_t(r)*nc+c]; if(!value) continue;
            out.write(reinterpret_cast<char*>(&r),4);out.write(reinterpret_cast<char*>(&c),4);out.write(reinterpret_cast<char*>(&value),8);nnz++;
        }
        out.close(); if(!out) throw std::runtime_error("Matrix write failed");
        std::ofstream audit(argv[7]);
        audit<<"{\"records\":"<<records<<",\"umi\":"<<total<<",\"observed_cells\":"<<observed.size()<<",\"nnz\":"<<nnz<<"}\n";
        if(!audit) throw std::runtime_error("Audit write failed");
    } catch(const std::exception& e) {std::cerr<<e.what()<<std::endl;return 1;}
}
