//
//  historyXmlWriter.h
//  VisIVOServer_New072012
//
//  Created by Fabio Vitello on 27/06/13.
//
//

#ifndef HISTORY_XML_WRITER_H
#define HISTORY_XML_WRITER_H


#include <string>
#include <vector>
#include <map>




class HistoryXmlWriter
{

    std::map<std::string, std::string> historyData;

public:
    HistoryXmlWriter(const char* histFile,const char* format,const char* out,const char* tableOrVolume,double comput[],double size[],
                                                        const char* login, const char* binaryHeader, float missing,float text, const char* endian,
                                                        const char* type, long unsigned int points,
                                                        const char* vo, const char* se, const char* lfnout, const char*  inputFile);
    
    HistoryXmlWriter(const char* histFile,const char* opName,std::map<std::string,std::string> appParameter, std::vector <std::string> outFilename);
    HistoryXmlWriter(const char* histFile, std::map<std::string,std::string> viewParameter,std::vector <std::string> outFilename);


    ~HistoryXmlWriter ( );
    
    void WriteXml(const char* histFile,std::map<std::string, std::string> h_data,std::vector <std::string> outFilename);
    
    void dump_map(const std::map<std::string, std::string>& map);


};


#endif

