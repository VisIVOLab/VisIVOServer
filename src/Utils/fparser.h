/***************************************************************************\
|* Function parser v2.83 by Warp                                           *|
|* -----------------------------                                           *|
|* Parses and evaluates the given function with the given variable values. *|
|* See fparser.txt for details.                                            *|
|*                                                                         *|
\***************************************************************************/

#ifndef ONCE_FPARSER_H_
#define ONCE_FPARSER_H_

#include <string>
#include <map>
#include <vector>

#ifdef FUNCTIONPARSER_SUPPORT_DEBUG_OUTPUT
#include <iostream>
#endif

class FunctionParser
{
public:
    enum ParseErrorType
    {
        SYNTAX_ERROR=0, MISM_PARENTH, MISSING_PARENTH, EMPTY_PARENTH,
        EXPECT_OPERATOR, OUT_OF_MEMORY, UNEXPECTED_ERROR, INVALID_VARS,
        ILL_PARAMS_AMOUNT, PREMATURE_EOS, EXPECT_PARENTH_FUNC,
        FP_NO_ERROR
    };


    int Parse(const std::string& Function, const std::string& Vars,
              bool useDegrees = false);
    const char* ErrorMsg() const;
    inline ParseErrorType GetParseErrorType() const { return parseErrorType; }

    double Eval(const double* Vars);
    inline int EvalError() const { return evalErrorType; }

    bool AddConstant(const std::string& name, double value);
    bool AddUnit(const std::string& name, double value);

    typedef double (*FunctionPtr)(const double*);

    bool AddFunction(const std::string& name,
                     FunctionPtr, unsigned paramsAmount);
    bool AddFunction(const std::string& name, FunctionParser&);
    void setNumberOfRow(unsigned int Rows){m_nOfRow=Rows;}
    void setNumberOfCol(unsigned int Cols){m_nOfCol=Cols;}
    void setArray(float **array,float* result){m_fArray=array;m_result=result;}

    void Optimize();


    FunctionParser();
    ~FunctionParser();

    // Copy constructor and assignment operator (implemented using the
    // copy-on-write technique for efficiency):
    FunctionParser(const FunctionParser&);
    FunctionParser& operator=(const FunctionParser&);


    void ForceDeepCopy();


#ifdef FUNCTIONPARSER_SUPPORT_DEBUG_OUTPUT
    // For debugging purposes only:
    void PrintByteCode(std::ostream& dest) const;
#endif



//========================================================================
private:
//========================================================================

// Private data:
// ------------
    ParseErrorType parseErrorType;
    int evalErrorType;

   float *m_result;
   int m_nOfRow;
   int m_nOfCol;
   float **m_fArray;

    struct Data
    {
        unsigned referenceCounter;

        int varAmount;
        bool useDegreeConversion, isOptimized;

        typedef std::map<std::string, unsigned> VarMap_t;
        VarMap_t Variables;

        typedef std::map<std::string, double> ConstMap_t;
        ConstMap_t Constants;
        ConstMap_t Units;

        VarMap_t FuncPtrNames;
        struct FuncPtrData
        {
            FunctionPtr ptr; unsigned params;
            FuncPtrData(FunctionPtr p, unsigned par): ptr(p), params(par) {}
        };
        std::vector<FuncPtrData> FuncPtrs;

        VarMap_t FuncParserNames;
        std::vector<FunctionParser*> FuncParsers;

        unsigned* ByteCode;
        unsigned ByteCodeSize;
        double* Immed;
        unsigned ImmedSize;
        double* Stack;
        unsigned StackSize;

        Data();
        ~Data();
        Data(const Data&);

        Data& operator=(const Data&); // not implemented on purpose
    };

    Data* data;
    unsigned evalRecursionLevel;

    // Temp data needed in Compile():
    unsigned StackPtr;
    std::vector<unsigned>* tempByteCode;
    std::vector<double>* tempImmed;


// Private methods:
// ---------------
    void CopyOnWrite();


    bool CheckRecursiveLinking(const FunctionParser*) const;

    bool ParseVars(const std::string& Vars,
                   std::map<std::string, unsigned>& dest);
    int VarNameType(const std::string&) const;
    Data::VarMap_t::const_iterator
    FindVariable(const char*, const Data::VarMap_t&) const;
    Data::ConstMap_t::const_iterator
    FindConstant(const char*, const Data::ConstMap_t&) const;
    int CheckForUnit(const char*, int) const;
    int CheckSyntax(const char*);
    int Compile(const char*);
    bool IsVariable(int);
    void AddCompiledByte(unsigned);
    void AddImmediate(double);
    void AddFunctionOpcode(unsigned);
    inline void incStackPtr();
    int CompileIf(const char*, int);
    int CompileFunctionParams(const char*, int, unsigned);
    int CompileElement(const char*, int);
    int CompilePossibleUnit(const char*, int);
    int CompilePow(const char*, int);
    int CompileUnaryMinus(const char*, int);
    int CompileMult(const char*, int);
    int CompileAddition(const char*, int);
    int CompileComparison(const char*, int);
    int CompileAnd(const char*, int);
    int CompileOr(const char*, int);
    int CompileExpression(const char*, int, bool=false);


    void MakeTree(void*) const;
};

#endif
