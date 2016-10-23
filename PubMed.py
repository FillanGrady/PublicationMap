from Bio import Entrez
import datetime
import os
Entrez.email = "fillangrady@gmail.com"
pickle_file_name = 'PickledLocations.txt'
output_file_name = 'Output.csv'
import googlemaps
gmaps = googlemaps.Client(key='AIzaSyA5Occ-Qhk8XEWV321gnGzAy3dWLHCawr0')
import geocoder
import time
import sys
import pickle
import pycountry


def month_to_int(month):
    import calendar
    try:
        return [k for k, v in enumerate(calendar.month_abbr) if v == month][0]
    except:
        raise ValueError("Month %s not found" % month)

class Rememberer:
    def __init__(self, f):
        self.f = f
        if os.path.exists(pickle_file_name):
            self.dict = pickle.load(open(pickle_file_name, "rb"))
            print "Loaded.  Size: %s" % (len(self.dict.keys()))
        else:
            self.dict = {}

    def execute(self, key):
        if key in self.dict:
            return self.dict[key]
        else:
            self.dict[key] = self.f(key)
            return self.execute(key)

    def pickle(self):
        print "Pickling.  Size: %s" % (len(self.dict.keys()))
        pickle.dump(self.dict, open(pickle_file_name, "wb"))


class Coordinate:
    def __init__(self, name):
        self.name = name
        self.coords = rememberer.execute(self.name)

    def __str__(self):
        return self.name + "," + str(self.coords[0]) + "," + str(self.coords[1])

    @staticmethod
    def geocode(name):
        global lookups
        lookups += 1
        g = geocoder.google(name)
        if g.latlng:
            global found_lookups
            found_lookups += 1
            return g.latlng
        else:
            return None, None


class Article:
    def __init__(self, record):
        self.title = record['MedlineCitation']['Article']['ArticleTitle']
        self.title = self.title.replace(",", "")
        self.locations = []
        for location in record['MedlineCitation']['Article']['AuthorList']:
            try:
                full_location = Article.remove_nonascii(location['AffiliationInfo'][0]['Affiliation'])
            except IndexError as e:
                continue
            loc = full_location.lower()
            loc = loc.replace("p.r. china", "china")
            loc = loc.replace("republic of korea", "korea")
            loc = loc.replace(";", ",")
            loc = loc.replace("-", ",")
            if "electronic" in loc:
                loc = loc[:loc.index("electronic")]
            loc = " ".join([x for x in loc.split(" ") if "@" not in x])
            institution = Article.get_institution(loc)
            country = Article.get_country(loc)
            loc = "%s %s" % (institution, country)
            loc = loc.strip()
            c = Coordinate(loc)
            if c.coords == (None, None):
                c = Coordinate(full_location)
                if c.coords == (None, None):
                    print "Not found: %s -> %s" % (full_location, loc)
            self.locations.append(c)

        r = record['MedlineCitation']['Article']
        try:
            date = record['MedlineCitation']['Article']['ArticleDate'][0]
            self.date = datetime.date(year=int(date['Year']), month=month_to_int(date['Month']), day=int(date['Day']))
        except IndexError as ie:
            self.date = None

    def __str__(self):
        return Article.remove_nonascii(self.title) + "," + \
               str(self.date) + "," + str(",".join([str(x) for x in self.locations]))

    @staticmethod
    def get_institution(full_location):
        if "university" in full_location:
            return [x.strip().strip(".") for x in full_location.split(",") if "university" in x][0]
        elif "institution" in full_location:
            return [x.strip().strip(".") for x in full_location.split(",") if "institution" in x][0]
        elif "institute" in full_location:
            return [x.strip().strip(".") for x in full_location.split(",") if "institute" in x][0]
        elif "hospital" in full_location:
            return [x.strip().strip(".") for x in full_location.split(",") if "hospital" in x][0]
        elif "center" in full_location:
            return [x.strip().strip(".") for x in full_location.split(",") if "center" in x][0]
        elif "laboratory" in full_location:
            return [x.strip().strip(".") for x in full_location.split(",") if "laboratory" in x][0]
        elif "college" in full_location:
            return [x.strip().strip(".") for x in full_location.split(",") if "college" in x][0]
        else:
            print "Couldn't determine institution: %s" % full_location
            return ""

    @staticmethod
    def get_country(full_location):
        words = full_location.split(",")
        for word in reversed(words):
            for country in countries:
                if country in word:
                    return word.strip().strip(".")
        return "usa"


    @staticmethod
    def remove_nonascii(s):
        return ''.join([i if ord(i) < 128 else '' for i in s])

rememberer = Rememberer(Coordinate.geocode)
countries = [x.name.lower() for x in pycountry.countries] + ["uk", "korea"]
country_abbreviations = [x.alpha3.lower() for x in pycountry.countries]


def get_links_term(term, retmax=2):
    links = Entrez.esearch(db="pubmed", retmax=retmax, term=term)
    record = Entrez.read(links)
    link_list = record[u'IdList']
    return link_list


def get_records(ids):
    handle = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
    records = Entrez.parse(handle=handle)
    return records

lookups = 0
found_lookups = 0
start = time.time()
topic = sys.argv[1]
retmax = sys.argv[2]
with open(output_file_name, "w") as f:
    g = get_records(get_links_term(topic, retmax=retmax))
    print "Time to get from pubmed: %s" % (str(time.time() - start))
    for i, a in enumerate(g):
        f.write(str(Article(a)) + os.linesep)
        if i % 50 == 0 and i > 0:
            rememberer.pickle()
print "Records: %s, Time: %s" % (retmax, str(time.time() - start))
print "Lookups: %s, Found Lookups: %s" % (str(lookups), str(found_lookups))