======================================
How to search for jobs
======================================

On the home page, navigate to the `Utils` section, located to the right of the logo.
Here, you'll find the `Job list` link.
Click it to be redirected to the job list page.

.. figure:: /_static/images/qeapp_home.png
   :align: center

The job list page displays a table of all jobs in the database, along with various filtering options to help you find specific jobs.
The table includes columns for the job ID, structure, creation time, state, label, relax type, and properties associated with each job.
Additionally, each row in the table provides links to `Delete` or `Inspect` the respective job via dynamically generated hyperlinks.

.. figure:: /_static/images/qeapp_job_search.png
   :align: center


Filters and Searches
--------------------
Users can filter the displayed data using various widgets:

   - **Properties Filter**: A set of checkboxes are generated for each unique property found across all jobs, allowing users to select jobs that match specific properties.
   - **Job State Dropdown**: Enables filtering based on the state of the job (e.g., finished, waiting, except, killed).
   - **Label Search Field**: A text input field where users can type a job label to narrow down the search results.
   - **Date Range Picker**: Two date pickers allow for filtering jobs between specific start and end dates.


The job list is updated dynamically based on the selected filters.

Delete job
-----------
In the delete page, the details of a specific job are shown.

.. figure:: /_static/images/qeapp_job_delete.png
   :align: center

1. **Display Job Information**:
   The system first displays the job's details, including its ID, type, label, description, and creation time.
   If there are any issues in loading the node details, an error message is displayed.

2. **Dependency Check**:
   Before deletion can proceed, the system checks for any dependent QEApp jobs linked to the node.
   If dependencies are found, the deletion is halted, and a warning is displayed.

3. **Confirmation Prompt**:
   If no dependencies are present, the system proceeds to ask for explicit user confirmation to ensure the action's seriousness and irreversibility is understood.
   The confirmation prompt is presented with a text input where the user must type "yes" to confirm.
   A button labeled "Delete Node" is provided to submit the confirmation.

By detailing each step and ensuring user confirmation, the system prevents accidental deletions of critical job data.
