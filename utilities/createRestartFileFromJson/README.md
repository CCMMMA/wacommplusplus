# WaComM++ utility for creating restart file starting from Json file

## Requirements (python3.6)



```bash
    python3.6 -m venv venv
    source venv/bin/activate
    pip3 install -r requirements.txt
```

## Usage

```bash
    python3 main.py <domain.nc> <features.json> <iDate>
```

### input

_domain.nc_: used to get lats and lons (WRF) \
_features.json_: stores info about area where create random particles \
_iDate_: YYYYMMDDZHH format

### output

NetCDF restart file containing particles

![img.png](img.png)

## Wacomm++

Inside configuration file (.json) set

```bash
...
    "restart": {
      "active": true,
      "restart_file": "restart_<iDate>.nc"
    },
    "sources": {
      "active": false,
      "sources_file": "sources.json"
    },
...
```