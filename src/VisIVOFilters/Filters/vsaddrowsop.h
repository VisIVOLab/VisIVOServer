#ifndef VSADDROWS_H
#define VSADDROWS_H

#include "vstableop.h"

class VSAddRowsOp : public VSTableOp
{
    std::vector< std::vector<float> > m_newRows;

public:
    VSAddRowsOp();
    ~VSAddRowsOp();

    void printHelp();
    bool execute();
};

#endif