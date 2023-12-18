# Filename: PathoPoleAnalyzer
FROM python:latest
COPY * ./
RUN python3 -m pip install -r requirements.txt

CMD [ "python"]
