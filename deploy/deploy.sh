#!/bin/bash
# Production Deployment Script for CriOS Dr. Crowe Coder
# Supports AWS, GCP, Azure, and DigitalOcean

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
REPO="ghcr.io/michaelcrowe11/crios"
VERSION=${VERSION:-latest}
ENVIRONMENT=${ENVIRONMENT:-production}

echo -e "${GREEN}═══════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}   CriOS Dr. Crowe Coder - Production Deployment${NC}"
echo -e "${GREEN}   194 PhD Agents Ready for Discovery${NC}"
echo -e "${GREEN}═══════════════════════════════════════════════════════${NC}\n"

# Check prerequisites
check_prerequisites() {
    echo -e "${YELLOW}Checking prerequisites...${NC}"
    
    # Check Docker
    if ! command -v docker &> /dev/null; then
        echo -e "${RED}Docker is not installed${NC}"
        exit 1
    fi
    
    # Check environment file
    if [ ! -f ".env.${ENVIRONMENT}" ]; then
        echo -e "${RED}.env.${ENVIRONMENT} file not found${NC}"
        exit 1
    fi
    
    echo -e "${GREEN}✓ Prerequisites met${NC}\n"
}

# Build and push Docker images
build_images() {
    echo -e "${YELLOW}Building Docker images...${NC}"
    
    # Backend
    docker build -t ${REPO}-backend:${VERSION} ./platform/backend
    docker push ${REPO}-backend:${VERSION}
    echo -e "${GREEN}✓ Backend image built and pushed${NC}"
    
    # Frontend
    docker build -t ${REPO}-frontend:${VERSION} ./platform/frontend
    docker push ${REPO}-frontend:${VERSION}
    echo -e "${GREEN}✓ Frontend image built and pushed${NC}\n"
}

# Deploy to AWS ECS
deploy_aws() {
    echo -e "${YELLOW}Deploying to AWS ECS...${NC}"
    
    # Update task definitions
    aws ecs register-task-definition --cli-input-json file://deploy/aws/task-definition.json
    
    # Update services
    aws ecs update-service \
        --cluster crios-cluster \
        --service crios-backend \
        --task-definition crios-backend:${VERSION} \
        --desired-count 4
    
    aws ecs update-service \
        --cluster crios-cluster \
        --service crios-frontend \
        --task-definition crios-frontend:${VERSION} \
        --desired-count 2
    
    echo -e "${GREEN}✓ Deployed to AWS ECS${NC}\n"
}

# Deploy to Google Cloud Run
deploy_gcp() {
    echo -e "${YELLOW}Deploying to Google Cloud Run...${NC}"
    
    # Deploy backend
    gcloud run deploy crios-backend \
        --image ${REPO}-backend:${VERSION} \
        --platform managed \
        --region us-central1 \
        --allow-unauthenticated \
        --min-instances 2 \
        --max-instances 10 \
        --memory 2Gi \
        --cpu 2 \
        --env-vars-file .env.${ENVIRONMENT}
    
    # Deploy frontend
    gcloud run deploy crios-frontend \
        --image ${REPO}-frontend:${VERSION} \
        --platform managed \
        --region us-central1 \
        --allow-unauthenticated \
        --min-instances 1 \
        --max-instances 5 \
        --memory 1Gi \
        --cpu 1
    
    echo -e "${GREEN}✓ Deployed to Google Cloud Run${NC}\n"
}

# Deploy to Azure Container Instances
deploy_azure() {
    echo -e "${YELLOW}Deploying to Azure Container Instances...${NC}"
    
    # Create resource group if not exists
    az group create --name crios-rg --location eastus
    
    # Deploy backend
    az container create \
        --resource-group crios-rg \
        --name crios-backend \
        --image ${REPO}-backend:${VERSION} \
        --cpu 2 \
        --memory 2 \
        --ports 8000 \
        --environment-variables-file .env.${ENVIRONMENT}
    
    # Deploy frontend
    az container create \
        --resource-group crios-rg \
        --name crios-frontend \
        --image ${REPO}-frontend:${VERSION} \
        --cpu 1 \
        --memory 1 \
        --ports 3000
    
    echo -e "${GREEN}✓ Deployed to Azure${NC}\n"
}

# Deploy with Docker Swarm
deploy_swarm() {
    echo -e "${YELLOW}Deploying with Docker Swarm...${NC}"
    
    # Initialize swarm if needed
    docker swarm init 2>/dev/null || true
    
    # Deploy stack
    docker stack deploy -c deploy/production.yml crios
    
    # Wait for services to be ready
    sleep 10
    
    # Show service status
    docker service ls | grep crios
    
    echo -e "${GREEN}✓ Deployed with Docker Swarm${NC}\n"
}

# Deploy to DigitalOcean App Platform
deploy_digitalocean() {
    echo -e "${YELLOW}Deploying to DigitalOcean...${NC}"
    
    doctl apps create --spec deploy/digitalocean/app-spec.yaml
    
    echo -e "${GREEN}✓ Deployed to DigitalOcean${NC}\n"
}

# Setup SSL with Let's Encrypt
setup_ssl() {
    echo -e "${YELLOW}Setting up SSL certificates...${NC}"
    
    # Install certbot if not present
    if ! command -v certbot &> /dev/null; then
        sudo apt-get update
        sudo apt-get install -y certbot python3-certbot-nginx
    fi
    
    # Get certificates
    sudo certbot --nginx -d crios.ai -d www.crios.ai -d api.crios.ai -d app.crios.ai \
        --non-interactive --agree-tos --email admin@crios.ai
    
    # Setup auto-renewal
    echo "0 0 * * * /usr/bin/certbot renew --quiet" | sudo crontab -
    
    echo -e "${GREEN}✓ SSL certificates configured${NC}\n"
}

# Run database migrations
run_migrations() {
    echo -e "${YELLOW}Running database migrations...${NC}"
    
    docker run --rm \
        --env-file .env.${ENVIRONMENT} \
        ${REPO}-backend:${VERSION} \
        python -m alembic upgrade head
    
    echo -e "${GREEN}✓ Database migrations complete${NC}\n"
}

# Setup monitoring
setup_monitoring() {
    echo -e "${YELLOW}Setting up monitoring...${NC}"
    
    # Deploy Prometheus
    docker service create \
        --name prometheus \
        --mount type=bind,source=$(pwd)/monitoring/prometheus.yml,target=/etc/prometheus/prometheus.yml \
        --publish 9090:9090 \
        prom/prometheus
    
    # Deploy Grafana
    docker service create \
        --name grafana \
        --env GF_SECURITY_ADMIN_PASSWORD=${GRAFANA_PASSWORD} \
        --publish 3001:3000 \
        grafana/grafana
    
    echo -e "${GREEN}✓ Monitoring setup complete${NC}\n"
}

# Health check
health_check() {
    echo -e "${YELLOW}Running health checks...${NC}"
    
    # Check backend
    if curl -f https://api.crios.ai/health > /dev/null 2>&1; then
        echo -e "${GREEN}✓ Backend is healthy${NC}"
    else
        echo -e "${RED}✗ Backend health check failed${NC}"
    fi
    
    # Check frontend
    if curl -f https://app.crios.ai > /dev/null 2>&1; then
        echo -e "${GREEN}✓ Frontend is healthy${NC}"
    else
        echo -e "${RED}✗ Frontend health check failed${NC}"
    fi
    
    echo ""
}

# Main deployment flow
main() {
    check_prerequisites
    
    # Parse deployment target
    case "${1}" in
        aws)
            build_images
            deploy_aws
            ;;
        gcp)
            build_images
            deploy_gcp
            ;;
        azure)
            build_images
            deploy_azure
            ;;
        digitalocean)
            build_images
            deploy_digitalocean
            ;;
        swarm)
            build_images
            deploy_swarm
            ;;
        *)
            echo -e "${YELLOW}Usage: ./deploy.sh [aws|gcp|azure|digitalocean|swarm]${NC}"
            echo -e "${YELLOW}Defaulting to Docker Swarm deployment...${NC}\n"
            build_images
            deploy_swarm
            ;;
    esac
    
    # Post-deployment tasks
    run_migrations
    setup_ssl
    setup_monitoring
    health_check
    
    echo -e "${GREEN}═══════════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}   Deployment Complete!${NC}"
    echo -e "${GREEN}   Platform: https://app.crios.ai${NC}"
    echo -e "${GREEN}   API: https://api.crios.ai${NC}"
    echo -e "${GREEN}   Docs: https://api.crios.ai/docs${NC}"
    echo -e "${GREEN}═══════════════════════════════════════════════════════${NC}\n"
}

# Run main function
main "$@"