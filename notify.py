import requests
import os
TG_TOKEN = os.environ['TG_TOKEN']
TG_CHAT_ID = os.environ['TG_CHAT_ID']
GROUP_NUMBER = os.environ['GROUP_NUMBER']
url = f"https://api.telegram.org/bot{TG_TOKEN}/sendMessage"
data = {
    "chat_id": TG_CHAT_ID,
    "text": f"Group: {GROUP_NUMBER} optimization completed"
}

response = requests.post(url, json=data)
print(response.json())