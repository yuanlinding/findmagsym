FROM python:3.11-slim

WORKDIR /app

# gcc/g++ needed for any deps that lack a pre-built wheel (spglib, spinspg)
RUN apt-get update \
    && apt-get install -y --no-install-recommends gcc g++ git \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY core.py .
COPY app.py .
COPY msg_list.cvs .

RUN adduser --disabled-password --gecos "" appuser
USER appuser

EXPOSE 8000
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000"]
