#!/bin/bash

# Prompt the user for yes/no input
read -p "Do you want to proceed? (yes/no): " response

# Convert the response to lowercase for case-insensitive comparison
response_lower=$(echo "$response" | tr '[:upper:]' '[:lower:]')

# Check the response and execute commands accordingly
if [[ "$response_lower" == "yes" ]]; then
    echo "Executing the command..."
    # Place your command here
    echo "Command executed."
elif [[ "$response_lower" == "no" ]]; then
    echo "Aborting."
else
    echo "Invalid response. Please enter 'yes' or 'no'."
fi
