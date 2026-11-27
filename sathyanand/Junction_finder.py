import requests

def get_junction_name(node_id):
    """Fetches location/junction name for a single OSM Node ID."""
    try:
        # 1. Fetch Node Lat/Lon from OSM API
        node_url = f"https://api.openstreetmap.org/api/0.6/node/{node_id}.json"
        res = requests.get(node_url, headers={'User-Agent': 'OSM-Reverse-Lookup'}).json()
        lat = res['elements'][0]['lat']
        lon = res['elements'][0]['lon']

        # 2. Reverse Geocode coordinates via Nominatim
        geo_url = f"https://nominatim.openstreetmap.org/reverse?format=json&lat={lat}&lon={lon}"
        geo_res = requests.get(geo_url, headers={'User-Agent': 'OSM-Reverse-Lookup'}).json()
        
        display_name = geo_res.get("display_name", "Unknown Location")
        return f"\n📍 Junction/Location Name:\n{display_name}\n(Coordinates: {lat}, {lon})"
    except Exception as e:
        return f"\nError resolving Node ID {node_id}: {e}"

def get_road_name(u_id, v_id):
    """Fetches road name connecting two OSM Node IDs (u and v)."""
    try:
        # Query Overpass API for way containing both nodes u and v
        overpass_url = "https://overpass-api.de/api/interpreter"
        query = f"[out:json];way(bn:{u_id})(bn:{v_id});out tags;"
        res = requests.post(overpass_url, data={'data': query}).json()
        elements = res.get('elements', [])

        if elements:
            tags = elements[0].get('tags', {})
            road_name = tags.get('name', 'Unnamed Road')
            highway_type = tags.get('highway', 'Unknown Type')
            return f"\n🛣️ Road Name: {road_name}\nHighway Type: {highway_type}"
        else:
            return f"\nNo direct OSM road way found connecting node {u_id} and {v_id}."
    except Exception as e:
        return f"\nError resolving edge ({u_id}, {v_id}): {e}"

if __name__ == "__main__":
    print("=== OpenStreetMap Reverse Lookup Tool ===")
    choice = input("Enter '1' to lookup a Junction (1 Node ID)\nEnter '2' to lookup a Road (2 Node IDs)\nChoice: ").strip()

    if choice == '1':
        node_id = input("\nEnter Junction Node ID: ").strip()
        print(get_junction_name(node_id))
    elif choice == '2':
        u_id = input("\nEnter Start Node ID (u): ").strip()
        v_id = input("Enter End Node ID (v): ").strip()
        print(get_road_name(u_id, v_id))
    else:
        print("Invalid selection. Please enter 1 or 2.")