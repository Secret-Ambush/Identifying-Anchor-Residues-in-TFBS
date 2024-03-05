from bs4 import BeautifulSoup
import json
import re

# Open the HTML file
with open('meme_out/meme.html', 'r') as file:
    html_content = file.read()

# Parse the HTML
soup = BeautifulSoup(html_content, 'html.parser')

script_tags = soup.find_all("script")
data_script = None
for script in script_tags:
    if 'var data' in script.text:
        data_script = script.text
        break

if data_script:
    json_str_match = re.search(r'var data = ({.*?});', data_script, re.DOTALL)
    if json_str_match:
        json_str = json_str_match.group(1)
        data = json.loads(json_str)
        pwm_section = data['motifs'][0]['pwm']
        print(pwm_section)
    else:
        print("JSON data not found in the script tag.")
else:
    print("Script tag containing 'var data' was not found.")
