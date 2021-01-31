//
// Created by Raffaele Montella on 08/12/20.
//

#include "Config.hpp"
#include "JulianDate.hpp"
#include <nlohmann/json.hpp>

// for convenience
using json = nlohmann::json;

Config::Config() {
    setDefault();
}

Config::Config(string &fileName): configFile(fileName) {
    setDefault();

    log4cplus::BasicConfigurator basicConfig;
    basicConfig.configure();
    logger = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT("WaComM"));

    LOG4CPLUS_DEBUG(logger, "Reading config file:" + fileName);

    if(fileName.substr(fileName.find_last_of(".") + 1) == "json") {
        // The configuration is a json
        loadFromJson(fileName);
    } else {
        // the configuration is a fortran style namelist
        loadFromNamelist(fileName);
    }

}

Config::~Config() {

}


void Config::setDefault() {
    // Is the model dry mode?
    dry = false;

    // Input time step in s
    _data.deltat=3600;

    // Integration time in s
    _data.dti=30;

    // To be clarified
    _data.tau0=86400.0;

    // Probability of particle surviving
    _data.survprob=1.0e-4;

    // Sedimentation velocity
    _data.sv=0;

    // Reduction Coefficient
    _data.crid=1;

    // File system root where the outputs will be saved ( default current directory )
    ncOutputRoot = "";

    timeStep = 3600;

    // Number of hours to calculate ( default 1)
    nHour = 1;

    // Starting time (default 0)
    startTime = 0;

    // Use restart (default false)
    useRestart = false;

    // Name of the file used for restarts
    restartFile = "WACOMM_rst_.txt";

    // Restart interval
    restartInterval = 3600;

    // Save history (to be used as restart) default true
    saveHistory = "none";

    // Default file name for the history
    historyRoot = "WACOMM_rst_";

    // Embedded history data into output (default false)
    embeddedHistory = false;

    // Use sources (default true)
    useSources = true;

    // Set the sources file name (defauly empty)
    sourcesFile = "";

    // Use random leap (defaukt true. use false for model testing versus other imlementations
    _data.random = true;

    // Save processed input files (default false)
    saveInput = false;

    // Root for saved input files.
    ncInputRoot="ocm_";
}

void Config::loadFromNamelist(const string &fileName) {
    setDefault();

    std::ifstream infile(fileName);
    std::string line;
    while (std::getline(infile, line))
    {
        // Trim the line
        line = Utils::trim(line," \t\r");

        // Select the section parser
        if (line == "&io") {
            // s starts with prefix &
            namelistParseIo(infile);
        } else if (line == "&chm") {
            // s starts with prefix &
            namelistParseChm(infile);
        } if (line == "&rst") {
            // s starts with prefix &
            namelistParseRst(infile);
        } if (line == "&hst") {
            // s starts with prefix &
            namelistParseHst(infile);
        }
    }
}

void Config::namelistParseIo(ifstream &infile) {
    std::string line;

    LOG4CPLUS_DEBUG(logger, "Parsing io section");

    // Start the input/output section
    while (std::getline(infile, line))
    {
        // Trim the line
        line = Utils::trim(line, " \t\r");

        // Check for the end of the section
        if (line=="/") { break; }

        vector<string> keyValues;
        Utils::tokenize(line,'=',keyValues);

        if (keyValues.size()==2) {
            string key=Utils::trim(keyValues.at(0)," \t\r");
            if (key == "nc_inputs") {
                // The list of the input files
                vector<string> ncInputs;
                Utils::tokenize(keyValues.at(1),',',ncInputs);

                for(string ncInput: ncInputs) {
                    ncInput=Utils::trim(ncInput," '");
                    this->ncInputs.push_back(ncInput);
                }
            } else if (key == "nc_output_root") {
                this->ncOutputRoot = Utils::trim(keyValues.at(1)," ',");
            } else if (key == "starttime") {
                this->startTime = stoi(keyValues.at(1));
            } else if (key == "nhour") {
                this->nHour = stoi(keyValues.at(1));
            } else if (key == "timestep") {
                this->timeStep = stod(keyValues.at(1));
            }
        }
    }
}

void Config::namelistParseChm(ifstream &infile) {
    std::string line;

    LOG4CPLUS_DEBUG(logger, "Parsing chm section");

    // Start the Chem section
    while (std::getline(infile, line))
    {
        // Trim the line
        line = Utils::trim(line, " \t\r");

        // Check for the end of the section
        if (line=="/") { break; }

        vector<string> keyValues;
        Utils::tokenize(line,'=',keyValues);

        if (keyValues.size()==2) {
            string key=Utils::trim(keyValues.at(0)," \t\r");
            if (key == "tau0") {
                this->_data.tau0 = stod(keyValues.at(1));
            } else if (key == "survprob") {
                this->_data.survprob = stod(keyValues.at(1));
            }
        }
    }
}

void Config::namelistParseRst(ifstream &infile) {
    std::string line;

    LOG4CPLUS_DEBUG(logger, "Parsing rst section");

    // Start the Restart section
    while (std::getline(infile, line))
    {
        // Trim the line
        line = Utils::trim(line, " \t\r");

        // Check for the end of the section
        if (line=="/") { break; }

        vector<string> keyValues;
        Utils::tokenize(line,'=',keyValues);

        if (keyValues.size()==2) {
            string key=Utils::trim(keyValues.at(0)," \t\r");
            if (key == "restart") {
                keyValues.at(1) = Utils::trim(keyValues.at(1), " \t");
                if (keyValues.at(1) == ".true.") {
                    this->useRestart=true;
                } else {
                    this->useRestart=false;
                }
            } else if (key == "restartfile") {
                this->restartFile = Utils::trim(keyValues.at(1)," '\t\r");;
            } else if (key == "interval") {
                this->restartInterval = stod(keyValues.at(1));
            }
        }
    }
}

void Config::namelistParseHst(ifstream &infile) {
    std::string line;

    LOG4CPLUS_DEBUG(logger, "Parsing hst section");

    // Start the History section
    while (std::getline(infile, line))
    {
        // Trim the line
        line = Utils::trim(line, " \t\r");

        // Check for the end of the section
        if (line=="/") { break; }

        vector<string> keyValues;
        Utils::tokenize(line,'=',keyValues);

        if (keyValues.size()==2) {
            string key=Utils::trim(keyValues.at(0)," \t\r");
            if (key == "history") {
                keyValues.at(1) = Utils::trim(keyValues.at(1), " \t");
                if (keyValues.at(1) == ".true.") {
                    this->saveHistory="text";
                }
            } else if (key == "historyfile") {
                this->historyRoot = Utils::trim(keyValues.at(1)," ',");
            }
        }
    }
}

double Config::Dti() const { return _data.dti; }

double Config::Deltat() const { return _data.deltat; }

double Config::Survprob() const {
    return _data.survprob;
}

double Config::Tau0() const {
    return _data.tau0;
}

double Config::SedimentationVelocity() const {
    return _data.sv;
}

bool Config::Dry() const {
    return dry;
}

void Config::Dry(bool value) {
    dry=value;
}

double Config::ReductionCoefficient() const {
    return _data.crid;
}

string Config::RestartFile() const {
    return restartFile;
}

bool Config::UseRestart() const {
    return useRestart;
}

string &Config::ConfigFile() {
    return configFile;
}

vector<string> &Config::NcInputs() {
    return ncInputs;
}

void Config::StartTimeIndex(int value) {
    startTime=value;
}

int Config::StartTimeIndex() {
    return startTime;
}

void Config::NumberOfInputs(int value) {
    nHour=value;
}

int Config::NumberOfInputs() {
    return nHour;
}



void Config::Random(bool value) {
    _data.random=value;
}

bool Config::Random() const {
    return _data.random;
}

config_data *Config::dataptr() {
    return &_data;
}

bool Config::UseSources() const {
    return useSources;
}

void Config::UseSources(bool value) {
    useSources = value;
}

void Config::UseRestart(bool value) {
    useRestart = value;
}

void Config::RestartFile(string value) {
    restartFile = value;
}

string Config::SourcesFile() const {
    return sourcesFile;
}

void Config::SourcesFile(string value) {
    sourcesFile=value;
}

string Config::NcOutputRoot() const {
    return ncOutputRoot;
}

void Config::NcOutputRoot(string value) {
    ncOutputRoot=value;
}

void  Config::SaveHistory(string value) {
    saveHistory=value;
}

void  Config::HistoryRoot(string value) {
    historyRoot = value;
}

string Config::SaveHistory() const {
    return saveHistory;
}

string Config::HistoryRoot() const {
    return historyRoot;
}

void Config::SaveInput(bool value) {
    saveInput = value;
}

bool Config::SaveInput() const {
    return saveInput;
}

string Config::NcInputRoot() const {
    return ncInputRoot;
}

void Config::NcInputRoot(string value) {
    ncInputRoot = value;
}

bool  Config::EmbeddedHistory() const {
    return embeddedHistory;
}

void Config::EmbeddedHistroy(bool value) {
    embeddedHistory = value;
}

string Config::OceanModel() const {
    return oceanModel;
}

void Config::OceanModel(string value) {
    oceanModel = value;
}

void Config::saveAsJson(const string &fileName) {

    Calendar calStart, calEnd;

    JulianDate::fromModJulian(julianStart, calStart);
    JulianDate::fromModJulian(julianEnd, calEnd);

    json simulation = {
            { "name", name },
            { "institution", institution },
            { "url", url},
            { "start", calStart.asNCEPdate()},
            { "end", calEnd.asNCEPdate() }
    };

    json io = {
            { "ocean_model", oceanModel},
            { "base_path", ncBasePath },
            { "nc_inputs", ncInputs },
            { "nc_output_root", ncOutputRoot },
            { "save_history", saveHistory },
            { "save_input", saveInput},
            { "nc_input_root", ncInputRoot}
    };

    json restart = {
            { "active", useRestart },
            { "restart_file", restartFile }
    };

    json sources = {
            { "active", useSources },
            { "sources_file", sourcesFile }
    };

    json physics = {
            { "tau0", _data.tau0 },
            { "survprob", _data.survprob },
            { "random", _data.random },
            { "sv", _data.sv },
            { "dti", _data.dti },
            { "deltat", _data.deltat },
            { "crid", _data.crid },
    };

    json config = {
            { "simulation", simulation},
            { "io", io},
            { "restart", restart},
            { "sources", sources},
            { "physics", physics},
    };

    // write prettified JSON to another file
    std::ofstream o(fileName);
    o << std::setw(4) << config << std::endl;
}

void Config::loadFromJson(const string &fileName) {
    setDefault();
    json config;
    std::ifstream i(fileName);
    i >> config;
    if (config.contains("simulation")) {
        json simulation=config["simulation"];
        if (simulation.contains("dry")) { dry = simulation["dry"]; }
        if (simulation.contains("name")) { name = simulation["name"]; }
        if (simulation.contains("institution")) { institution = simulation["institution"]; }
        if (simulation.contains("url")) { url = simulation["url"]; }
        if (simulation.contains("start")) {
            Calendar calStart(simulation["start"]);
            julianStart=JulianDate::toModJulian(calStart);
        }
        if (simulation.contains("end")) {
            Calendar calEnd(simulation["end"]);
            julianEnd=JulianDate::toModJulian(calEnd);
        }
    }
    if (config.contains("io")) {
        json io=config["io"];
        if (io.contains("embedded_history")) { embeddedHistory = io["embedded_history"]; }
        if (io.contains("save_history")) { saveHistory = io["save_history"]; }
        if (io.contains("history_root")) {  historyRoot= io["history_root"]; }
        if (io.contains("save_input")) { saveInput = io["save_input"]; }
        if (io.contains("nc_input_root")) { ncInputRoot = io["nc_input_root"]; }
        if (io.contains("base_path")) { ncBasePath = io["base_path"]; }
        if (io.contains("ocean_model")) { oceanModel = io["ocean_model"]; }
        if (io.contains("nc_output_root")) { ncOutputRoot = io["nc_output_root"]; }
        if (io.contains("timestep")) { timeStep = io["timestep"]; }
        if (io.contains("nc_inputs") && io["nc_inputs"].is_array()) {
            for (auto ncInput:io["nc_inputs"]) {
                string file=ncInput;
                this->ncInputs.push_back(ncBasePath+"/"+file);
            }
        }
    }
    if (config.contains("restart")) {
        json restart=config["restart"];
        if (restart.contains("active")) { useRestart = restart["active"]; }
        if (restart.contains("restart_file")) { restartFile = restart["restart_file"]; }
    }
    if (config.contains("sources")) {
        json sources=config["sources"];
        if (sources.contains("active")) { useSources = sources["active"]; }
        if (sources.contains("sources_file")) { sourcesFile = sources["sources_file"]; }
    }
    if (config.contains("physics")) {
        json physics=config["physics"];
        if (physics.contains("tau0")) { _data.tau0 = physics["tau0"]; }
        if (physics.contains("crid")) { _data.crid = physics["crid"]; }
        if (physics.contains("deltat")) { _data.deltat = physics["deltat"]; }
        if (physics.contains("dti")) { _data.dti = physics["dti"]; }
        if (physics.contains("sv")) { _data.sv = physics["sv"]; }
        if (physics.contains("random")) { _data.random = physics["random"]; }
        if (physics.contains("survprob")) { _data.survprob = physics["survprob"]; }
    }
}

