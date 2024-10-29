FROM python:3.12-slim

# install package
ENV DEBIAN_FRONTEND=noninteractive
RUN \
    sed -i 's/deb.debian.org/mirrors.tuna.tsinghua.edu.cn/' /etc/apt/sources.list.d/debian.sources && \
    sed -i 's/security.debian.org/mirrors.tuna.tsinghua.edu.cn/' /etc/apt/sources.list.d/debian.sources && \
    sed -i 's/security-cdn.debian.org/mirrors.tuna.tsinghua.edu.cn/' /etc/apt/sources.list.d/debian.sources && \
    apt update && \
    apt install --no-install-recommends -y \
    unzip zip default-jre locales \
    r-base r-base-dev && \
    apt clean && rm -rf /var/lib/apt/lists/*

# set chinese fonts
RUN sed -ie 's/# zh_CN.UTF-8 UTF-8/zh_CN.UTF-8 UTF-8/g' /etc/locale.gen
RUN locale-gen
ENV LANG=zh_CN.UTF-8
ADD tools/fonts.tar.xz /usr/share/fonts/
RUN fc-cache -vf
RUN fc-list -v | grep zh

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PYTHONIOENCODING=UTF-8
RUN pip3 config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple/ && \
    pip3 config set install.trusted-host pypi.tuna.tsinghua.edu.cn

# ADD https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz /usr/local/src

# ENV PATH="$PATH:/usr/local/src/ncbi-blast-2.15.0+/bin"
# RUN mkdir -p /usr/local/src/blast && \
#     tar xvf /usr/local/src/ncbi-blast-2.15.0+-x64-linux.tar.gz -C /usr/local/src/ && \
#     rm -f /usr/local/src/ncbi-blast-2.15.0+-x64-linux.tar.gz

ADD . /app/

WORKDIR /app

RUN \
    --mount=type=cache,mode=0755,target=/root/.cache/pip \
    python -m pip install --upgrade pip && \
    pip install -r requirements.txt

# 通过清华镜像源安装
RUN R -e "install.packages('foreach', repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"
RUN R -e "install.packages('doParallel', repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"

# ENTRYPOINT ["python3", "app.py"]

RUN python down_ref.py
CMD python run_xaem.py
