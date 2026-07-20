Main Work Flow Steps:
If MasterTable.csv, UserTable.csv, and ExpertxUser.csv exists:
1. Load those three tables from GitHub
2. Load ExpertxUserOnward
3. Connect Google API (Steps at End)
4. Run each block accordingly
5. Push new UserTable.csv and ExpertxUser.csv to GitHub
6. Move used student csvs out of "Please upload the file containing your photometry results ("YBphotometry_results_instructorID.csv") here. (File responses)" 

WARNING! WARNING!
DO NOT RUN MORE THAN ONCE WITH SAME DATASETS OR THE DATA WILL BE SKEWED

AverageMaskTable:
1. Connect Google API (Steps At End)
2. First block pulls files from the folder with all the Expert CSVs and stors in a local folder that code pulls from
4. The next block loads an individual file from drive that will be used for header formatting
3. The next block is then ran which averages the photometry and masks which keeps any pixel that 50% of the expert masks agree on. This can be changed by changing mask_quality variable. 
4. This is all saved to a csv (MasterTable.csv) with the new masks, photometry values, standard error, and a count of photometry values used.
5. All tables are created within the current working directory

CheckingqualInitial:
1. Connect Google API (Steps At End)
2. Second block loads data from MasterTable and UsedTables (A list of all used files)
3. Next block pulls files from the folder with all the student values and saves to a local folder to pull from. 
4. Next block defines functions that will be used in block 6
5. Block 5 takes master table and creates a Boolean mask from each wavelength of each yellowball. It is then stored and saved to ExpertBool. This steamlines the process so each mask doesn't need to create a Boolean mask each time the code runs.
6. Next block quality checks the students mask by photometry and mask quality. If the photometry isn't within 60% of the expert photometry it isn't accepted. If the student mask does match the expert maks more than 60% the photom value isn't included. These values can be changed in mask_check and photom_check.
7. The code then averages in the accepted student photometry values and increases the count accordingly. There is also a count kept of that of photometry values excluded. The standard error then is updated for the new average.
8. This is all saved to a csv (ExpertxUser.csv) with the new masks, photometry values, standard error, and a count of photometry values used.
9. The last block then updates the UsedTables.csv to include all tables used 
10. All tables are created within the current working directory

CheckingqualOnward:
1. Connect Google API (Steps At End)
2. Second block loads data from MasterTable, UsedTables (A list of all used files), and ExpertxUser
3. Next block pulls files from the folder with all the new student values and saves to a local folder to pull from. 
4. Next block defines functions that will be used in block 5
5. Next block quality checks the students mask by photometry and mask quality. If the photometry isn't within 60% of the expert photometry it isn't accepted. If the student mask does match the expert maks more than 60% the photom value isn't included. These values can be changed in mask_check and photom_check.
6. The code then averages in the accepted student photometry values to the ExpertxUser values and increases the count accordingly. There is also a count kept of that of photometry values excluded. The standard error then is updated for the new average. 
7. This is all saved to a csv (ExpertxUser.csv) with the new masks, photometry values, standard error, and a count of photometry values used.
8. The last block then updates the UsedTables.csv to include all tables used 
9. All tables are created within the current working directory

Google API Steps: 
1. Go To https://console.cloud.google.com/
2. Sign in with the google account with the files
3. Click select a project
4. Click new project
5. Name it and click create
6. Then go to APIs and Services
7. Click library
8. Search for google drive api and click it 
9. Click enable
10. Go back to apis and services then to Oauth consent screen
11. Click external the click create
12. Put in your project name then your email twice
13. Then go to audience under the oauth consent page
14. Click add users and then add your gmail address
15. Then enter api and services then go to credentials
16. Click create credentials then oauth client ID
17. Click desktop app and name it then click create
18. Download the JSON file and move it to project folder
19. Install packages:
a. pip install google-api-python-client
b. pip install google-auth-httplib2
c. pip install google-auth-oauthlib
20. Replace current json file in "flow" with new one’s name

