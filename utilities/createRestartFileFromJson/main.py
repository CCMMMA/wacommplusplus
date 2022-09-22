import sys
import json
import random
import netCDF4
import numpy as np
from datetime import datetime
from shapely.geometry import Polygon, Point
from tqdm import tqdm


# Program to generate restart file starting from bboxes stored in json file


# Pass from lat, lon and depth to indices (i, j, k)
def deplatlon2kji(depth, lat, lon):
    # print((np.abs(np.asarray(latDomain) - lat)).argmin())
    k = (np.abs(np.asarray(depthDomain) - depth)).argmin()
    j = (np.abs(np.asarray(latDomain) - lat)).argmin()
    i = (np.abs(np.asarray(lonDomain) - lon)).argmin()

    return k, j, i


# Create NetCDF file
def createNetCDF():
    print("Creating NetCDF ...")

    ncdstfile = netCDF4.Dataset("restart_" + iDate + ".nc", "w", format="NETCDF4")
    ncdstfile.createDimension("particles", size=count)
    ncdstfile.createDimension("particle_time", size=particle_time)

    timeVar = ncdstfile.createVariable("particle_time", "i4", "particle_time")
    timeVar.description = "Time since initialization"
    timeVar.long_name = "time since initialization"
    timeVar.units = "seconds since 1968-05-23 00:00:00"
    timeVar.calendar = "gregorian"
    timeVar.field = "time, scalar, series"
    timeVar[:] = timestamp

    idVar = ncdstfile.createVariable("id", "f4", "particles")
    idVar.long_name = "id"
    idVar.cf_role = "trajectory_id"
    idVar[:] = id

    depthVar = ncdstfile.createVariable("depth", "f4", ("particles", "particle_time"))
    depthVar.description = "depth"
    depthVar.long_name = "depth"
    depthVar.units = "meters"
    depthVar[:] = deps

    lonVar = ncdstfile.createVariable("longitude", "f4", ("particles", "particle_time"))
    lonVar.description = "Longitude"
    lonVar.long_name = "longitude"
    lonVar.units = "degrees_east"
    lonVar[:] = lons

    latVar = ncdstfile.createVariable("latitude", "f4", ("particles", "particle_time"))
    latVar.description = "Latitude"
    lonVar.long_name = "latitude"
    latVar.units = "degrees_north"
    latVar[:] = lats

    kVar = ncdstfile.createVariable("k", "f4", ("particles", "particle_time"))
    kVar.long_name = "depth fractional index"
    kVar[:] = K

    JVar = ncdstfile.createVariable("j", "f4", ("particles", "particle_time"))
    JVar.long_name = "latitude fractional index"
    JVar[:] = J

    IVar = ncdstfile.createVariable("i", "f4", ("particles", "particle_time"))
    IVar.long_name = "longitude fractional index"
    IVar[:] = I

    HealthVar = ncdstfile.createVariable("health", "f4", ("particles", "particle_time"))
    HealthVar.long_name = "health of the particle"
    HealthVar.units = "1e-9"
    HealthVar[:] = health

    HealthVar = ncdstfile.createVariable("age", "f4", ("particles", "particle_time"))
    HealthVar.long_name = "age of the particle"
    HealthVar.units = "seconds since emission"
    HealthVar[:] = age

    timeVar = ncdstfile.createVariable("time", "i4", ("particles", "particle_time"))
    timeVar.description = "emission time of the particle"
    timeVar.long_name = "time since initialization"
    timeVar.units = "seconds since 1968-05-23 00:00:00"
    timeVar.calendar = "gregorian"
    timeVar[:] = time

    ncdstfile.close()


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 " + sys.argv[0] + " domain.nc features.geojson iDate")
        exit(0)

    # Parsing arguments
    domainFile = netCDF4.Dataset(sys.argv[1])
    featuresFile = open(sys.argv[2])
    iDate = sys.argv[3]

    # Load features from json file
    features = json.load(featuresFile)["features"]

    latDomain = domainFile.variables['latitude']
    lonDomain = domainFile.variables['longitude']
    depthDomain = domainFile.variables['depth']

    count = 0
    particle_time = 1

    id = []
    lats = []
    lons = []
    deps = []
    J = []
    I = []
    K = []
    health = []
    age = []
    time = []

    # Calculate timestamp since 1968-05-23 00:00:00
    actualDate = datetime.strptime(iDate, "%Y%m%dZ%H%M")
    date1968 = datetime(1968, 5, 23)
    timestamp = (actualDate - date1968).total_seconds()

    # For each feature
    for feature in features:
        # Get coordinates
        coordinates = feature["geometry"]["coordinates"][0]
        particles = feature["properties"]["particles"]

        # Create polygon from coordinates
        X, Y = zip(*coordinates)
        polygon = Polygon(zip(X, Y))
        minx, miny, maxx, maxy = polygon.bounds

        print(f'Generating random particles ({particles}) ...')

        # For each particle to generate
        for idx in tqdm(range(0, particles)):
            # Generate random position (lat and lon) inside polygon
            inside = False
            while not inside:
                lon = random.uniform(minx, maxx)
                lat = random.uniform(miny, maxy)
                dep = feature["properties"]["depth"]

                particle = Point(lon, lat)
                # Check if generated particle is inside area
                if polygon.contains(particle):
                    # Get k, j and i from dep, lat, lon
                    k, j, i = deplatlon2kji(dep, lat, lon)

                    id.append(count)
                    lats.append(lat)
                    lons.append(lon)
                    deps.append(dep)
                    J.append(j)
                    I.append(i)
                    K.append(k)
                    health.append(1)
                    age.append(0)
                    time.append(timestamp)

                    count += 1
                    inside = True

    createNetCDF()
    print("End!")
