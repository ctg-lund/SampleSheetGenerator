# CTG Samplesheet Generator
## Introduction
This application is used for samplesheet generation and parsing at the next generation sequencing core facility Center for Translational Genomics (CTG). It will, given a csv file, create a samplesheet that can both be used for demultilexing with bcl-convert, and further downstream analysis with our meta pipelines [Yggdrasil](https://github.com/ctg-lund/Yggdrasil).
## Quick start
In order to start the application, make sure that you have docker, docker compose, and git installed. Then run the following: 
```
git clone https://github.com/ctg-lund/SampleSheetGenerator
cd SampleSheetGenerator
docker-compose up
```
Then navigate to localhost:8082/samplesheet.

## Technical description
The tech stack is based on three main technologies:
* Flask
  * For building the actual web application which can process user data and return a validated samplesheet in pure python (with a little jinja sprinkled on top)
* Gunicorn
  * For a production grade wsgi, instead of flask's built-in wsgi which isn't meant to be in production
* Pandas
  * For easier parsing and handling of csv data, used in the "backend" of the application.

The application is built with the intention to be deployed behind a proxy apache2 server. For more information on how to deploy on a apache server, [look here](docs/deployment.md)!

# Updates
Form updated
main script takes index kit from form
- add indexes to samplesheet
- add fields to replace a well with the index pair of another well
- strip unnecessary info from illumina v2
- start from a simple csv table with just sample names filled in