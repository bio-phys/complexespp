//////////////////////////////////////////////////////
/// MIT License, 2022
/// Author: Bérenger Bramas
/// See https://gitlab.inria.fr/bramas/clsimple
//////////////////////////////////////////////////////

#ifndef CLSIMPLE_HPP
#define CLSIMPLE_HPP

#include <string>
#include <vector>
#include <memory>
#include <optional>
#include <regex>
#include <sstream>
#include <typeinfo>
#include <type_traits>
#include <limits>
#include <cassert>
#include <set>

class CLsimple{
    template <class ParamType>
    static ParamType Convert(const std::string& str, bool* outFlag = nullptr){
        std::istringstream iss(str, std::istringstream::in);
        ParamType value;
        iss >> value;
        if(outFlag){
            (*outFlag) = bool(iss.eof());
        }
        return value;
    }

    template <class ParamType>
    static std::string TypeToStr(){
        if constexpr(std::is_same<ParamType,bool>::value){
            return "Boolean";
        }
        if constexpr(std::numeric_limits<ParamType>::is_integer){
            return "Integer";
        }
        if constexpr(std::is_floating_point<ParamType>::value){
            return "Real number";
        }
        if constexpr(std::is_same<ParamType,std::string>::value){
            return "String";
        }
        return typeid(std::remove_reference_t<ParamType>).name();
    }


    enum ParamTypes{
        NoArg,
        Multi,
        Single
    };

    class AbstractParam{
        const std::vector<std::string> _keys;
        const std::string _help;
        const ParamTypes _paramType;
        const int _mandatoryGroup;

    public:
        AbstractParam(std::vector<std::string> inKeys,
                      std::string inHelp,
                      const ParamTypes inType,
                      const int inMandatoryGroup):
            _keys(std::move(inKeys)), _help(std::move(inHelp)),
            _paramType(inType), _mandatoryGroup(inMandatoryGroup){
            assert(_keys.size());
        }

        virtual ~AbstractParam(){}

        std::vector<std::string> getKeys() const{
            return _keys;
        }

        std::string getHelp() const{
            return _help;
        }

        bool isSingle() const {
            return _paramType == Single;
        }

        bool isMulti() const {
            return _paramType == Multi;
        }

        bool isNoArg() const {
            return _paramType == NoArg;
        }

        bool isMandatory() const{
            return _mandatoryGroup != NotMandatory;
        }

        int mandatoryGroup() const{
            return _mandatoryGroup;
        }

        virtual bool applyValue(const std::string& inValue) = 0;
        virtual void applyDefault() = 0;
        virtual std::string getTypeStr() const = 0;
    };

    template <class ParamType>
    class MultiParam : public AbstractParam{
        std::optional<std::reference_wrapper<std::vector<ParamType>>> _variable;
        const std::vector<ParamType> _default;

    public:
        MultiParam(std::vector<std::string> inKeys,
                   std::string inHelp,
                   const int inMandatoryGroup,
                   std::optional<std::reference_wrapper<std::vector<ParamType>>> inVariable,
                   std::vector<ParamType> inDefault):
            AbstractParam(std::move(inKeys), std::move(inHelp), Multi, inMandatoryGroup),
            _variable(std::move(inVariable)),
            _default(std::move(inDefault)){}

        bool applyValue(const std::string& inValue) final{
            if constexpr (std::is_same<ParamType, bool>::value){
                bool convertionOkInt;
                auto valueInt = Convert<int>(inValue, &convertionOkInt);
                if(_variable && convertionOkInt){
                    _variable->get().push_back(bool(valueInt));
                    return true;
                }
                bool convertionOkStr;
                auto valueStr = Convert<std::string>(inValue, &convertionOkStr);
                if(_variable && convertionOkStr){
                    _variable->get().push_back(valueStr == "TRUE" || valueStr == "true" || valueStr == "True");
                    return true;
                }
                return false;
            }
            else{
                bool convertionOk;
                auto value = Convert<ParamType>(inValue, &convertionOk);
                if(_variable && convertionOk){
                    _variable->get().push_back(value);
                    return true;
                }
                return false;
            }
        }
        void applyDefault() final{
            if(_variable){
                _variable->get() = _default;
            }
        }
        std::string getTypeStr() const final{
            return "List of " + TypeToStr<ParamType>() + "s";
        }
    };


    template <class ParamType>
    class Param : public AbstractParam{
        std::optional<std::reference_wrapper<ParamType>> _variable;
        const ParamType _default;

    public:
        Param(std::vector<std::string> inKeys,
              std::string inHelp,
              const int inMandatoryGroup,
              std::optional<std::reference_wrapper<ParamType>> inVariable,
              ParamType inDefault):
            AbstractParam(std::move(inKeys), std::move(inHelp), Single, inMandatoryGroup),
            _variable(std::move(inVariable)),
            _default(std::move(inDefault)){}

        bool applyValue(const std::string& inValue) final{
            if constexpr (std::is_same<ParamType, bool>::value){
                bool convertionOkInt;
                auto valueInt = Convert<int>(inValue, &convertionOkInt);
                if(_variable && convertionOkInt){
                    _variable->get() = bool(valueInt);
                    return true;
                }
                bool convertionOkStr;
                auto valueStr = Convert<std::string>(inValue, &convertionOkStr);
                if(_variable && convertionOkStr){
                    _variable->get() = (valueStr == "TRUE" || valueStr == "true" || valueStr == "True");
                    return true;
                }
                return false;
            }
            else{
                bool convertionOk;
                auto value = Convert<ParamType>(inValue, &convertionOk);
                if(_variable && convertionOk){
                    _variable->get() = value;
                    return true;
                }
                _variable->get() = _default;
                return false;
            }
        }
        void applyDefault() final{
            if(_variable){
                _variable->get() = _default;
            }
        }
        std::string getTypeStr() const final{
            return TypeToStr<ParamType>();
        }
    };

    class ParamNoArg : public AbstractParam{
    public:
        ParamNoArg(std::vector<std::string> inKeys,
              std::string inHelp,
                   const int inMandatoryGroup):
            AbstractParam(std::move(inKeys), std::move(inHelp), NoArg, inMandatoryGroup){}

        bool applyValue(const std::string& inValue) final{
            return true;
        }
        void applyDefault() final{
        }
        std::string getTypeStr() const final{
            return "(No argument)";
        }
    };

    static bool StrSeemsAKey(const std::string& str){
        std::regex keyFormat("-{1,2}(\\D|[0-9]\\D)");
        return std::regex_search(str, keyFormat);
    }

    static std::string StrToKey(const std::string& str){
        std::regex keyValueFormat("-{1,2}([^=]+).*");
        std::smatch match;
        std::regex_search(str, match, keyValueFormat);
        return match[1];
    }

    static bool IsKeyValueFormat(const std::string& str){
        std::regex keyValueFormat("--[^=]+=.+");
        return std::regex_search(str, keyValueFormat);
    }

    static std::pair<std::string,std::string> SplitKeyValue(const std::string& str){
        std::regex keyValueFormat("--([^=]+)=(.+)");
        std::smatch match;
        std::regex_search(str, match, keyValueFormat);
        const std::string key = match[1];
        const std::string value = match[2];
        return std::pair<std::string,std::string>(std::move(key), std::move(value));
    }

    template <class ParamType>
    void processParam(ParamType&& param, bool* parseIsOK = nullptr, std::set<int>* usedFields = nullptr) {
        const auto& keys = param->getKeys();
        int pos = -1;
        for(const auto& key : keys){
            const int testPos = getKeyPos(key);
            if(testPos != -1){
                pos = testPos;
                break;
            }
        }

        if(pos == -1){
            param->applyDefault();
        }
        else{
            if(usedFields){
                (*usedFields).insert(pos);
            }
            if(param->isNoArg()){
                // Do nothing
            }
            else if(IsKeyValueFormat(_argv[pos])){
                const auto keyValue = SplitKeyValue(_argv[pos]);
                const bool res = param->applyValue(std::get<1>(keyValue));
                if(!res){
                    (*errorStr) << "[ERROR] Error in parsing the value for arg \""
                             << _argv[pos] << "\"\n";
                }
                if(parseIsOK){
                    (*parseIsOK) &= res;
                }
            }
            else if(pos+1 != int(_argv.size())){
                if(param->isMulti()){
                    int idxVal = pos + 1;
                    while(idxVal != int(_argv.size()) && !StrSeemsAKey(_argv[idxVal])){
                        const bool res = param->applyValue(_argv[idxVal]);
                        if(!res){
                            (*errorStr) << "[ERROR] Error in parsing the value \"" << _argv[idxVal] << "\" for arg \""
                                     << _argv[pos] << "\"\n";
                        }
                        if(parseIsOK){
                            (*parseIsOK) &= res;
                        }
                        if(usedFields){
                            (*usedFields).insert(idxVal);
                        }
                        idxVal += 1;
                    }
                    if(idxVal == pos + 1){
                        param->applyDefault();

                        (*errorStr) << "[ERROR] Cannot get value, but we expect one, for arg \""
                                 << _argv[pos] << "\"\n";

                        if(parseIsOK){
                            (*parseIsOK) = false;
                        }
                    }
                }
                else{
                    const bool res = param->applyValue(_argv[pos+1]);
                    if(!res){
                        (*errorStr) << "[ERROR] Error in parsing the value \"" << _argv[pos+1] << "\" for arg \""
                                 << _argv[pos] << "\"\n";
                    }
                    if(parseIsOK){
                        (*parseIsOK) &= res;
                    }
                    if(usedFields){
                        (*usedFields).insert(pos+1);
                    }
                }
            }
            else{
                param->applyDefault();
                (*errorStr) << "[ERROR] Cannot get value, but we expect one, for arg \""
                         << _argv[pos] << "\"\n";
                if(parseIsOK){
                    (*parseIsOK) = false;
                }
            }
        }
    }

    const std::string _title;
    std::string _exec;

    std::vector<std::string> _argv;
    std::shared_ptr<std::vector<std::unique_ptr<AbstractParam>>> _params;

    bool _failsIfInvalid;
    bool _acceptUnregisteredParams;
    bool _isValid;

    std::shared_ptr<std::ostringstream> errorStr;

public:
    static const int NotMandatory = -1;

    CLsimple(const std::string inTitle,
             const int argc, const char *const argv[],
             const bool inFailsIfInvalid = true,
             const bool inAcceptUnregisteredParams = false)
        : _title(std::move(inTitle)),
          _params(new std::vector<std::unique_ptr<AbstractParam>>),
          _failsIfInvalid(inFailsIfInvalid),
          _acceptUnregisteredParams(inAcceptUnregisteredParams),
          _isValid(true),
          errorStr(new std::ostringstream){
        if(argc){
            _exec = argv[0];
            _argv.reserve(argc-1);
            for(int idxArg = 1 ; idxArg < argc ; ++idxArg){
                _argv.emplace_back(argv[idxArg]);
            }
        }
    }

    CLsimple(const CLsimple&) = default;
    CLsimple& operator=(const CLsimple&) = default;

    bool failsIfInvalid() const{
        return _failsIfInvalid;
    }

    bool acceptUnregisteredParams() const{
        return _acceptUnregisteredParams;
    }

    bool isValid() const{
        return _isValid;
    }

    int getKeyPos(const std::string& inKey) const{
        for(int idxArg = 0 ; idxArg < int(_argv.size()) ; ++idxArg){
            if(StrSeemsAKey(_argv[idxArg]) && StrToKey(_argv[idxArg]) == inKey){
                return idxArg;
            }
        }
        return -1;
    }

    bool hasKey(const std::string& inKey) const{
        return getKeyPos(inKey) != -1;
    }

    bool hasOneOfKeys(const std::vector<std::string>& inKeys) const{
        for(const auto& key : inKeys){
            const int pos = getKeyPos(key);
            if(pos != -1){
                return true;
            }
        }
        return false;
    }

    bool hasOneOfKeys(std::initializer_list<std::string> inKeys) const{
        std::vector<std::string> keys = inKeys;
        return hasOneOfKeys(keys);
    }

    template <class ParamType>
    std::vector<ParamType> getValues(std::vector<std::string> inKeys,
                  std::vector<ParamType> inDefaultValue = std::vector<ParamType>(),
                  bool* outOk = nullptr) const{
        std::vector<ParamType> variable;
        std::optional<std::reference_wrapper<std::vector<ParamType>>> variableRef(variable);
        MultiParam<ParamType> newParam(std::move(inKeys), "", false,
                                                    std::move(variableRef),
                                                    std::move(inDefaultValue));
        processParam(&newParam, outOk);
        return variable;
    }

    template <class ParamType>
    ParamType getValue(std::vector<std::string> inKeys,
                  ParamType inDefaultValue = ParamType(),
                  bool* outOk = nullptr) const{
        ParamType variable;
        std::optional<std::reference_wrapper<ParamType>> variableRef(variable);
        Param<ParamType> newParam(std::move(inKeys), "", false,
                                                    std::move(variableRef),
                                                    std::move(inDefaultValue));
        processParam(&newParam, outOk);
        return variable;
    }

    bool parse(){
        bool parseIsOK = true;
        std::set<int> usedFields;

        // Process param values
        for(auto& param : *_params){
            processParam(param, &parseIsOK, &usedFields);
        }
        if(_failsIfInvalid){
            _isValid = parseIsOK;
            if(!_acceptUnregisteredParams && usedFields.size() != _argv.size()){
                (*errorStr) << "[ERROR] The following arguments are invalid:\n";
                for(int idx = 0 ; idx < int(_argv.size()) ; ++idx){
                    if(usedFields.find(idx) == usedFields.end()){
                        (*errorStr) << "[ERROR] \"" << _argv[idx] << "\"\n";
                    }
                }
                _isValid = false;
            }
        }
        // Ensure that one mandatory group is used
        std::unordered_map<int, bool> groupStates;
        for(auto& param : *_params){
            if(param->isMandatory()){
                if(groupStates.find(param->mandatoryGroup()) == groupStates.end()){
                    groupStates[param->mandatoryGroup()] = true;
                }
                const auto& keys = param->getKeys();
                groupStates[param->mandatoryGroup()] &= hasOneOfKeys(keys);
            }
        }
        if(groupStates.size()){
            bool oneGroupIsOk = false;
            for(auto iter : groupStates){
                oneGroupIsOk = iter.second;
                if(oneGroupIsOk){
                    break;
                }
            }
            if(!oneGroupIsOk){
                (*errorStr) << "[ERROR] Some arguments are mandatory, but none of them were provided...\n";
            }

            _isValid &= oneGroupIsOk;
        }
        return _isValid;
    }

    template <class ParamType>
    void addMultiParameter(std::vector<std::string> inKeys, std::string inHelp,
                           std::optional<std::reference_wrapper<std::vector<ParamType>>> inVariable = std::nullopt,
                           std::vector<ParamType> inDefaultValue = std::vector<ParamType>(),
                           const int inMandatoryGroup = NotMandatory){
        std::unique_ptr<AbstractParam> newParam(new MultiParam<ParamType>(
                                                    std::move(inKeys),
                                                    std::move(inHelp),
                                                    inMandatoryGroup,
                                                    std::move(inVariable),
                                                    std::move(inDefaultValue)
                                                    ));
        (*_params).emplace_back(std::move(newParam));
    }

    template <class ParamType>
    void addParameter(std::vector<std::string> inKeys, std::string inHelp,
                      std::optional<std::reference_wrapper<ParamType>> inVariable = std::nullopt,
                      ParamType inDefaultValue = ParamType(),
                      const int inMandatoryGroup = NotMandatory){
        std::unique_ptr<AbstractParam> newParam(new Param<ParamType>(
                                                    std::move(inKeys),
                                                    std::move(inHelp),
                                                    inMandatoryGroup,
                                                    std::move(inVariable),
                                                    std::move(inDefaultValue)
                                                    ));
        (*_params).emplace_back(std::move(newParam));
    }

    void addParameterNoArg(std::vector<std::string> inKeys, std::string inHelp,
                           const int inMandatoryGroup = NotMandatory){
        std::unique_ptr<AbstractParam> newParam(new ParamNoArg(
                                                    std::move(inKeys),
                                                    std::move(inHelp),
                                                    inMandatoryGroup
                                                    ));
        (*_params).emplace_back(std::move(newParam));
    }

    template <class StreamClass>
    void printHelp(StreamClass& inStream, const bool inPrintError = true) const{
        auto join = [](auto&& elements, auto&& delimiter) -> std::string {
            std::string res;
            for(const auto& elm : elements){
                if(std::size(res)){
                    res += delimiter;
                }
                res += elm;
            }
            return res;
        };

        inStream << "[HELP] \n";
        inStream << "[HELP] " << _title << "\n";
        inStream << "[HELP] \n";
        if(!_isValid){
            inStream << "[HELP] \n";
            inStream << "[HELP] Currant parameters are not valid...\n";
            inStream << "[HELP] \n";
        }
        for(auto& param : *_params){
            inStream << "[HELP] Parameter names: {" << join(param->getKeys(),", ") << "}\n";
            inStream << "[HELP]  - Description: " << param->getHelp() << "\n";
            inStream << "[HELP]  - Type: " << param->getTypeStr() << "\n";
            if(param->isMandatory()){
                inStream << "[HELP]  - Is mandatory (group " << param->mandatoryGroup() << ")\n";
            }
            inStream << "[HELP]\n";
        }
        inStream << "[HELP]" << std::endl;


        inStream << "[HELP] Command line with all args: " << _exec << " ";
        for(auto& param : *_params){
            assert(param->getKeys().size());
            if(param->isMulti()){
                inStream << "--" << param->getKeys()[0] << " [" << param->getTypeStr() << "] ";
            }
            else if(param->isSingle()){
                inStream << "--" << param->getKeys()[0] << "=[" << param->getTypeStr() << "] ";
            }
            else{
                inStream << "--" << param->getKeys()[0] << " ";
            }
        }
        inStream << std::endl;
        inStream << "[HELP]" << std::endl;

        if(inPrintError){
            printError(inStream);
        }
    }

    template <class StreamClass>
    void printError(StreamClass& inStream) const{
        const std::string errors = (*errorStr).str();
        if(errors.length()){
            inStream << "[ERROR] There are some errors...\n";
            inStream << errors << std::endl;;
        }
    }

    template <class ValType>
    static ValType GetMapping(const std::string& inKey, const std::vector<std::pair<std::string,ValType>>& inMapping,
                              const ValType defaultVal = ValType()){
        for(auto& mp : inMapping){
            if(mp.first == inKey){
                return mp.second;
            }
        }
        return defaultVal;
    }

    template <class ValType>
    static ValType GetMapping(const std::string& inKey,
                              std::initializer_list<std::pair<std::string, ValType>> inMapping,
                              const ValType defaultVal = ValType()){
        for(auto& mp : inMapping){
            if(mp.first == inKey){
                return mp.second;
            }
        }
        return defaultVal;
    }
};


#endif
