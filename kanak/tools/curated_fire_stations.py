# Author: Kanak Varma
"""
Curated Official Fire Stations for Bengaluru (KSFES) and Chennai (TNFRS).
Sourced directly from:
  - Karnataka State Fire & Emergency Services (KSFES) GIS/KML directory via OpenCity.in & KSRSAC
  - Greater Chennai Corporation (GCC) & Tamil Nadu Fire & Rescue Services (TNFRS) official directory
"""

CHENNAI_FIRE_STATIONS = [
    # ── Chennai North Division ──
    {"name": "Esplanade Fire Station", "lat": 13.086478, "lon": 80.286064, "address": "Esplanade, Rajaji Salai, Chennai"},
    {"name": "Korukkupet Fire Station", "lat": 13.111800, "lon": 80.272100, "address": "Kathbada, Korukkupet, Chennai"},
    {"name": "Royapuram Fire Station", "lat": 13.109200, "lon": 80.294600, "address": "East Kalmandapam Road, Royapuram, Chennai"},
    {"name": "Secretariat Fire Station", "lat": 13.079800, "lon": 80.287100, "address": "Fort St. George Complex, Chennai"},
    {"name": "Sembium Fire Station", "lat": 13.124900, "lon": 80.239100, "address": "Madurai Sami Madam St, Sembiyam, Perambur, Chennai"},
    {"name": "Tondiarpet Fire Station", "lat": 13.124500, "lon": 80.288200, "address": "Thiruvottiyur High Road, Tondiarpet, Chennai"},
    {"name": "V.O.C. Nagar Fire Station", "lat": 13.116500, "lon": 80.264200, "address": "V.O.C. Nagar, Tondiarpet, Chennai"},
    {"name": "Vyasarpadi Fire Station", "lat": 13.118400, "lon": 80.259400, "address": "Sathyamoorthy Nagar Main Rd, Vyasarpadi, Chennai"},
    {"name": "Washermanpet Fire Station", "lat": 13.102440, "lon": 80.274401, "address": "Maniyakaran St, Old Washermanpet, Chennai"},
    {"name": "Moolakothalam Fire Station", "lat": 13.105508, "lon": 80.278214, "address": "Moolakothalam, Old Washermanpet, Chennai"},

    # ── Chennai South Division ──
    {"name": "Ashok Nagar Fire Station", "lat": 13.033600, "lon": 80.210400, "address": "Ashok Pillar Main Rd, Zaffarkhanpet, Chennai"},
    {"name": "Guindy Fire Station", "lat": 13.009400, "lon": 80.212500, "address": "SIDCO Industrial Estate, Alandur Rd, Guindy, Chennai"},
    {"name": "Mylapore Fire Station", "lat": 13.035800, "lon": 80.267300, "address": "Kutchery Road, Mylapore, Chennai"},
    {"name": "Rajbhavan Fire Station", "lat": 13.011200, "lon": 80.222800, "address": "Sardar Patel Road, Raj Bhavan, Guindy, Chennai"},
    {"name": "Saidapet Fire Station", "lat": 13.021100, "lon": 80.222500, "address": "Anna Salai, Little Mount / Saidapet, Chennai"},
    {"name": "Teynampet Fire Station", "lat": 13.047447, "lon": 80.249119, "address": "Anna Salai, Teynampet, Chennai"},
    {"name": "T. Nagar Fire Station", "lat": 13.038500, "lon": 80.231200, "address": "South Usman Road / Bus Terminus, T. Nagar, Chennai"},
    {"name": "Thiruvanmiyur Fire Station", "lat": 12.983100, "lon": 80.259400, "address": "East Coast Road (ECR), Thiruvanmiyur, Chennai"},
    {"name": "Velachery Fire Station", "lat": 12.977200, "lon": 80.222100, "address": "Tiruvalluvar St, Annai Indira Nagar, Velachery, Chennai"},
    {"name": "Tambaram Fire Station", "lat": 12.925400, "lon": 80.117200, "address": "GST Road, State Training Centre, Tambaram, Chennai"},
    {"name": "Tambaram Sanatorium Fire Station", "lat": 12.943636, "lon": 80.133596, "address": "GST Road, Tambaram Sanatorium, Chennai"},

    # ── Chennai Central Division ──
    {"name": "Kilpauk Fire Station", "lat": 13.078200, "lon": 80.242500, "address": "Kilpauk Garden Road, Kilpauk, Chennai"},
    {"name": "Vepery Fire Station", "lat": 13.083500, "lon": 80.264700, "address": "EVK Sampath Salai, Vepery, Chennai"},
    {"name": "Egmore Fire Station (Directorate)", "lat": 13.068633, "lon": 80.261333, "address": "Pantheon Road, Egmore, Chennai"},
    {"name": "Triplicane Fire Station", "lat": 13.064500, "lon": 80.278500, "address": "Wallajah Road / Anna Salai, Triplicane, Chennai"},
    {"name": "Koyambedu Fire Station", "lat": 13.068900, "lon": 80.191800, "address": "Market A Road, Wholesale Market Complex, Koyambedu, Chennai"},
    {"name": "J.J. Nagar Fire Station", "lat": 13.084100, "lon": 80.176200, "address": "Valayapathi Salai, Mogappair East, Chennai"},
    {"name": "Villivakkam Fire Station", "lat": 13.107200, "lon": 80.205400, "address": "4th Main Road, SIDCO Nagar, Villivakkam, Chennai"},

    # ── Chennai Suburban & Industrial Corridor ──
    {"name": "Ambattur Fire Station", "lat": 13.102100, "lon": 80.158200, "address": "3rd Main Rd, Ambattur Industrial Estate, Chennai"},
    {"name": "Avadi Fire Station", "lat": 13.111394, "lon": 80.106389, "address": "TNHB, Avadi, Chennai"},
    {"name": "Madhavaram Fire Station", "lat": 13.155000, "lon": 80.241200, "address": "Milk Colony Road, Madhavaram, Chennai"},
    {"name": "Manali Fire Station", "lat": 13.163600, "lon": 80.258600, "address": "CPCL Industrial Belt, Manali, Chennai"},
    {"name": "Thiruvottiyur Fire Station", "lat": 13.160200, "lon": 80.301500, "address": "Rajakadai, T.H. Road, Thiruvottiyur, Chennai"},
    {"name": "Ennore Fire Station", "lat": 13.208500, "lon": 80.321800, "address": "Kathivakkam High Road, Ennore, Chennai"},
    {"name": "Poonamallee Fire Station", "lat": 13.049200, "lon": 80.094100, "address": "Trunk Road, Poonamallee, Chennai"},
    {"name": "Redhills Fire Station", "lat": 13.190200, "lon": 80.187400, "address": "Sothupakkam Rd / GNT Rd, Red Hills, Chennai"},
    {"name": "Sholinganallur Fire Station", "lat": 12.901500, "lon": 80.228500, "address": "Rajiv Gandhi Salai (OMR), Sholinganallur, Chennai"},
    {"name": "Siruseri SIPCOT Fire Station", "lat": 12.830800, "lon": 80.222300, "address": "SIPCOT IT Park, OMR, Siruseri, Chennai"},
    {"name": "Chennai Airport Fire Station", "lat": 12.992503, "lon": 80.166022, "address": "Airports Authority of India, Meenambakkam, Chennai"},
]

BENGALURU_FIRE_STATIONS = [
    # ── Central / BBMP Core ──
    {"name": "Bengaluru Secretariat Fire Station", "lat": 12.977766, "lon": 77.589087, "address": "Vidhana Soudha, Bengaluru"},
    {"name": "High Grounds Fire Station", "lat": 12.980723, "lon": 77.580304, "address": "Millers Road, High Grounds, Bengaluru"},
    {"name": "Ulsoor / Sivanchetti Fire Station", "lat": 12.981200, "lon": 77.618500, "address": "Annaswamy Mudaliar Road, Halasuru, Bengaluru"},
    {"name": "Bengaluru South Fire Station (HQ)", "lat": 12.972409, "lon": 77.610860, "address": "Field Marshal Cariappa Rd, Ashok Nagar, Bengaluru"},
    {"name": "Bengaluru North Fire Station", "lat": 12.962526, "lon": 77.563082, "address": "Sirsi Circle, Mysore Road, Chamarajpet, Bengaluru"},
    {"name": "Rajajinagar Fire Station", "lat": 12.980473, "lon": 77.548617, "address": "1st Main Rd, Industrial Town, Rajajinagar, Bengaluru"},

    # ── North & North-East Radial ──
    {"name": "Hebbal Fire Station", "lat": 13.048455, "lon": 77.618139, "address": "Bellary Road, Hebbal, Bengaluru"},
    {"name": "Thanisandra Fire Station", "lat": 13.051200, "lon": 77.625400, "address": "Manyata Tech Park / Thanisandra Main Rd, Bengaluru"},
    {"name": "Yelahanka Fire Station", "lat": 13.125133, "lon": 77.600948, "address": "BBMP Complex, Yelahanka, Bengaluru"},
    {"name": "Banaswadi Fire Station", "lat": 13.011564, "lon": 77.645920, "address": "Dodda Banaswadi Main Rd, Bengaluru"},
    {"name": "ITI Fire Station", "lat": 13.003798, "lon": 77.678640, "address": "Old Madras Road, Dooravani Nagar, Bengaluru"},
    {"name": "Hoskote Fire Station", "lat": 13.066050, "lon": 77.779000, "address": "National Highway 75, Hoskote, Bengaluru"},

    # ── East & IT Corridor ──
    {"name": "Mahadevapura Fire Station", "lat": 12.998590, "lon": 77.702669, "address": "Outer Ring Road, Mahadevapura, Bengaluru"},
    {"name": "Whitefield Fire Station", "lat": 12.987356, "lon": 77.732640, "address": "EPIP Zone, Whitefield, Bengaluru"},
    {"name": "Sarjapur Road Fire Station", "lat": 12.916863, "lon": 77.673745, "address": "Sarjapur Main Road, Bellandur / Kaikondrahalli, Bengaluru"},
    {"name": "Electronics City Fire Station", "lat": 12.835533, "lon": 77.662695, "address": "Phase 1, Electronic City, Hosur Road, Bengaluru"},

    # ── South & South-West ──
    {"name": "Jayanagara Fire Station", "lat": 12.935181, "lon": 77.600854, "address": "Bannerghatta Road / 9th Block Jayanagar, Bengaluru"},
    {"name": "Banashankari Fire Station", "lat": 12.922952, "lon": 77.561995, "address": "100 Feet Ring Road, Banashankari, Bengaluru"},
    {"name": "Hulimavu Fire Station", "lat": 12.872264, "lon": 77.598791, "address": "Bannerghatta Main Road, Hulimavu, Bengaluru"},
    {"name": "Anjanapura Fire Station", "lat": 12.859067, "lon": 77.563267, "address": "Anjanapura 80 Feet Road, Bengaluru"},

    # ── West & North-West ──
    {"name": "Nagarbhavi Fire Station", "lat": 12.935527, "lon": 77.509925, "address": "Bangalore University Campus, Mysore Rd, Nagarbhavi, Bengaluru"},
    {"name": "Sunkadakatte Fire Station", "lat": 12.974212, "lon": 77.421651, "address": "Magadi Main Road, Sunkadakatte, Bengaluru"},
    {"name": "Peenya Fire Station", "lat": 13.028053, "lon": 77.520669, "address": "Peenya 1st Stage, Industrial Area, Bengaluru"},
    {"name": "Yeshwantpur Fire Station", "lat": 13.022558, "lon": 77.549758, "address": "Tumkur Road, Yeshwantpur, Bengaluru"},
    {"name": "Nelamangala Fire Station", "lat": 13.099800, "lon": 77.389300, "address": "Tumkur Road Bypass, Nelamangala, Bengaluru"},
]

