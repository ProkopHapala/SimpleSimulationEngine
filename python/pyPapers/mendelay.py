from mendeley import Mendeley
from mendeley.session import MendeleySession

# Replace with your client_id and client_secret
client_id = 'your_client_id'
client_secret = 'your_client_secret'
redirect_uri = 'http://localhost:8080'  # This should match the redirect URI you registered

# Initialize the Mendeley object
mendeley = Mendeley(client_id, client_secret, redirect_uri=redirect_uri)

# Get the authorization URL
auth = mendeley.start_authorization_code_flow()
print(f"Go to the following URL to authorize the application: {auth.get_login_url()}")

# Copy the authorization code from the URL and paste it here
authorization_code = input("Enter the authorization code: ")

# Authenticate using the authorization code
session = auth.authenticate(authorization_code)

# You can now use the session to interact with the Mendeley API

# Example: List all documents in your Mendeley library
documents = session.library.all()

for document in documents.iter():
    print(f"{document.title} by {document.authors}")

# Example: Get folders
folders = session.folders.iter()

for folder in folders:
    print(f"Folder: {folder.name}")




# # List Documents:
# documents = session.library.all()
# for document in documents.iter():
#     print(f"{document.title} by {document.authors}")

# #Get Folders:
# folders = session.folders.iter()
# for folder in folders:
#     print(f"Folder: {folder.name}")

# #Get Documents in a Folder:
# folder_id = 'your_folder_id'  # Replace with the actual folder ID
# folder_documents = session.folders.get(folder_id).documents()
# for document in folder_documents.iter():
#     print(f"{document.title} by {document.authors}")

# #Search for a Document:
# search_query = 'machine learning'
# search_results = session.catalog.search(search_query).list()
# for doc in search_results:
#     print(f"{doc.title} by {doc.authors}")