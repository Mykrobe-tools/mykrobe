FROM python:3.6
RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app
RUN pip install --upgrade pip
COPY . /usr/src/app
RUN pip install -r requirements.txt 
RUN python setup.py install
CMD mykrobe --help