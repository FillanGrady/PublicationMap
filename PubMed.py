from Bio import Entrez
import datetime
import os
Entrez.email = "fillangrady@gmail.com"
pickle_file_name = 'PickledLocations.txt'
output_file_name = 'Output.csv'
key1 = 'AIzaSyA5Occ-Qhk8XEWV321gnGzAy3dWLHCawr0'
key2 = 'AIzaSyC69qprzE_k_WggYOdHQSddnfJgcASIrS8'
key3 = 'AIzaSyDQzC6_3HLdCrv55J9CDZYR0t7OeFhTqK4'
key4 = 'AIzaSyBoiY3tbFEkqYA2Rbr8p4zteBIC_IyYKdk'
import time
import sys
import pickle
import pycountry
from googleplaces import GooglePlaces, GooglePlacesError


def month_to_int(month):
    import calendar
    try:
        return [k for k, v in enumerate(calendar.month_abbr) if v == month][0]
    except:
        raise ValueError("Month %s not found" % month)


class GooglePlaces_MultipleKeys:
    def __init__(self, key_list):
        self.key_list = key_list
        self.index = 0
        self.google_place = GooglePlaces(key_list[self.index])

    def text_search(self, name):
        try:
            global lookups
            lookups += 1
            query_result = self.google_place.text_search(name)
            location = None
            for place in query_result.places:
                location = place.geo_location
            if location is None:
                global found_lookups
                found_lookups += 1
                return 0, 0
            else:
                return float(location['lat']), float(location['lng'])
        except GooglePlacesError as ge:
            self.index += 1
            if self.index > len(self.key_list) - 1:
                raise GooglePlacesError("Out of keys")
            print "Switching Keys"
            self.google_place = GooglePlaces(self.key_list[self.index])
            return self.text_search(name)

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
    def __init__(self, full_name):
        self.full_name = full_name
        self.remove_emails()
        self.wo_department = self.create_wo_department()
        self.institution, self.country = self.create_institution_country()
        self.institution_country = "%s %s" % (self.institution, self.country)
        self.city_country = "%s %s" % self.create_city_country()
        self.coords = rememberer.execute(self.institution_country)
        if self.coords == (0, 0):
            self.coords = rememberer.execute(self.city_country)
            if self.coords == (0, 0):
                self.coords = rememberer.execute(self.wo_department)
                if self.coords == (0, 0):
                    self.coords = rememberer.execute(self.full_name)
                    if self.coords == (0, 0):
                        print "-----------------------"
                        print self.full_name
                        print self.wo_department
                        print self.institution_country
                        print "________________________"

    def remove_emails(self):
        self.full_name = self.full_name.lower()
        self.full_name = self.full_name.replace(";", ",")
        self.full_name = self.full_name.replace("-", ",")
        if "electronic" in self.full_name:
            self.full_name = self.full_name[:self.full_name.index("electronic")]
        self.full_name = " ".join([x for x in self.full_name.split(" ") if "@" not in x])

    def create_wo_department(self):
        return ",".join([x for x in self.full_name.split(",") if not
            any(dpt in x for dpt in ["department"])])

    def create_institution_country(self):
        institution = Coordinate.get_institution(self.full_name)
        country = Coordinate.get_country(self.full_name)
        return institution, country

    def create_city_country(self):
        city = self.full_name.split(",")[-2]
        country = Coordinate.get_country(self.full_name)
        return city, country

    def __str__(self):
        return "%s,%s,%.3f,%.3f" % (self.full_name.replace(",", ""), self.country, self.coords[0], self.coords[1])

    @staticmethod
    def geocode(name):
        if name == "":
            return 0, 0
        return google_places.text_search(name)

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
            return ""

    @staticmethod
    def get_country(full_location):
        words = full_location.split(",")
        for word in reversed(words):
            for country in countries:
                if country in word:
                    return word.strip().strip(".")
        return ""


class Article:
    def __init__(self, record):
        self.title = Article.remove_nonascii(record['MedlineCitation']['Article']['ArticleTitle'])
        self.title = self.title.replace(",", "")
        self.locations = []
        try:
            for location in record['MedlineCitation']['Article']['AuthorList']:
                try:
                    full_location = Article.remove_nonascii(location['AffiliationInfo'][0]['Affiliation'])
                except (IndexError, KeyError):
                    continue
                c = Coordinate(full_location)
                self.locations.append(c)
        except KeyError as ke:
            pass
        try:
            date = record['MedlineCitation']['Article']['ArticleDate'][0]
            self.date = datetime.date(year=int(date['Year']), month=month_to_int(date['Month']), day=int(date['Day']))
        except IndexError as ie:
            self.date = None

    def __str__(self):
        return "%s,%s,%s" % (self.title, str(self.date), ",".join([str(x) for x in self.locations]))

    @staticmethod
    def remove_nonascii(s):
        return ''.join([i if ord(i) < 128 else '' for i in s])

rememberer = Rememberer(Coordinate.geocode)
countries = [x.name.lower() for x in pycountry.countries] + ["uk", "korea", "usa"]
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
google_places = GooglePlaces_MultipleKeys([key1, key2, key3, key4])
with open(output_file_name, "w") as f:
    g = get_records(get_links_term(topic, retmax=retmax))
    print "Time to get from pubmed: %s" % (str(time.time() - start))
    for i, a in enumerate(g):
        f.write(str(Article(a)) + os.linesep)
        if i % 50 == 0 and i > 0:
            rememberer.pickle()
print "Records: %s, Time: %s" % (retmax, str(time.time() - start))
print "Lookups: %s, Found Lookups: %s" % (str(lookups), str(found_lookups))
rememberer.pickle()