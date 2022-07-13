library(googlesheets4)
# Set authentication token to be stored in a folder called `.secrets`
options(gargle_oauth_cache = ".secrets")

# Authenticate manually
gs4_auth()


# If successful, the previous step stores a token file.
# Check that a file has been created with:
list.files(".secrets/")

# Check that the non-interactive authentication works by first deauthorizing:
gs4_deauth()

# Authenticate using token. If no browser opens, the authentication works.
gs4_auth(cache = ".secrets", email = "k.pham@yale.edu")


ss <- gs4_get("https://docs.google.com/spreadsheets/d/1cHCdfy1WiIOFFJDEpo3U97KaiPep0Tv_dvtIEEcnE3c/edit#gid=93956359")
sheet_append(ss, data.frame(time=Sys.time()))
