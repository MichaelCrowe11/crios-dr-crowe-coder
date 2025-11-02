#!/bin/bash

# CriOS Dr. Crowe Coder - Railway Deployment Script
# This script automates the deployment of all services to Railway

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if Railway CLI is installed
check_railway_cli() {
    log_info "Checking Railway CLI installation..."
    if ! command -v railway &> /dev/null; then
        log_error "Railway CLI is not installed!"
        echo "Install via: npm i -g @railway/cli"
        exit 1
    fi
    log_success "Railway CLI is installed"
}

# Check if logged in to Railway
check_railway_login() {
    log_info "Checking Railway authentication..."
    if ! railway whoami &> /dev/null; then
        log_warning "Not logged in to Railway"
        log_info "Running: railway login"
        railway login
    fi
    log_success "Authenticated with Railway"
}

# Create or link Railway project
setup_railway_project() {
    log_info "Setting up Railway project..."

    if [ ! -f ".railway/config.json" ]; then
        log_warning "No Railway project found"
        echo -n "Would you like to create a new project or link an existing one? (create/link): "
        read choice

        if [ "$choice" == "create" ]; then
            railway init
        elif [ "$choice" == "link" ]; then
            railway link
        else
            log_error "Invalid choice"
            exit 1
        fi
    fi

    log_success "Railway project configured"
}

# Add PostgreSQL
add_postgres() {
    log_info "Checking PostgreSQL database..."

    if railway variables | grep -q "POSTGRES_HOST"; then
        log_success "PostgreSQL already provisioned"
    else
        log_info "Adding PostgreSQL database..."
        railway add --plugin postgresql
        log_success "PostgreSQL added successfully"
    fi
}

# Add Redis
add_redis() {
    log_info "Checking Redis cache..."

    if railway variables | grep -q "REDIS_HOST"; then
        log_success "Redis already provisioned"
    else
        log_info "Adding Redis cache..."
        railway add --plugin redis
        log_success "Redis added successfully"
    fi
}

# Set environment variables
set_environment_variables() {
    log_info "Setting up environment variables..."

    # Check if .env.railway exists
    if [ ! -f ".env.railway" ]; then
        log_warning ".env.railway not found"
        log_info "Please create .env.railway based on .env.railway.template"
        echo -n "Continue without setting custom variables? (y/n): "
        read continue

        if [ "$continue" != "y" ]; then
            exit 0
        fi
    else
        log_info "Loading environment variables from .env.railway..."
        # Railway CLI will automatically load from .env.railway
    fi

    # Prompt for critical variables if not set
    if ! railway variables | grep -q "ANTHROPIC_API_KEY"; then
        log_warning "ANTHROPIC_API_KEY not set"
        echo -n "Enter your Anthropic API key (or press Enter to skip): "
        read api_key

        if [ ! -z "$api_key" ]; then
            railway variables set ANTHROPIC_API_KEY="$api_key"
            log_success "ANTHROPIC_API_KEY set"
        fi
    fi

    # Set JWT_SECRET if not set
    if ! railway variables | grep -q "JWT_SECRET"; then
        log_info "Generating JWT_SECRET..."
        jwt_secret=$(openssl rand -base64 32)
        railway variables set JWT_SECRET="$jwt_secret"
        log_success "JWT_SECRET generated and set"
    fi
}

# Deploy services
deploy_services() {
    log_info "Starting deployment of all services..."

    # Backend API
    log_info "Deploying Backend API..."
    railway up --service crios-backend || log_warning "Backend deployment may have issues - check logs"

    # Frontend
    log_info "Deploying Frontend..."
    railway up --service crios-frontend || log_warning "Frontend deployment may have issues - check logs"

    # Celery Worker
    log_info "Deploying Celery Worker..."
    railway up --service crios-worker || log_warning "Worker deployment may have issues - check logs"

    # Core API (optional)
    log_info "Deploying Core API..."
    railway up --service crios-core-api || log_warning "Core API deployment may have issues - check logs"

    log_success "All services deployment initiated"
}

# Health check
run_health_checks() {
    log_info "Running health checks (waiting 60s for services to start)..."
    sleep 60

    # Get service URLs
    backend_url=$(railway status --service crios-backend --json | jq -r '.deployments[0].url' 2>/dev/null || echo "")
    frontend_url=$(railway status --service crios-frontend --json | jq -r '.deployments[0].url' 2>/dev/null || echo "")

    if [ ! -z "$backend_url" ]; then
        log_info "Checking backend health at $backend_url/health"
        if curl -f -s "$backend_url/health" > /dev/null; then
            log_success "Backend is healthy"
        else
            log_warning "Backend health check failed"
        fi
    fi

    if [ ! -z "$frontend_url" ]; then
        log_info "Checking frontend at $frontend_url"
        if curl -f -s "$frontend_url" > /dev/null; then
            log_success "Frontend is accessible"
        else
            log_warning "Frontend health check failed"
        fi
    fi
}

# Display deployment summary
show_summary() {
    echo ""
    echo "=========================================="
    log_success "Deployment Complete!"
    echo "=========================================="
    echo ""
    log_info "Access your services:"
    echo ""
    railway status
    echo ""
    log_info "View logs:"
    echo "  railway logs --service crios-backend --follow"
    echo ""
    log_info "Manage environment variables:"
    echo "  railway variables"
    echo ""
    log_info "For more information, see RAILWAY_DEPLOYMENT.md"
    echo ""
}

# Main execution
main() {
    echo ""
    echo "=========================================="
    echo "  CriOS Dr. Crowe Coder - Railway Deploy  "
    echo "=========================================="
    echo ""

    check_railway_cli
    check_railway_login
    setup_railway_project
    add_postgres
    add_redis
    set_environment_variables
    deploy_services
    run_health_checks
    show_summary

    log_success "Deployment script completed!"
}

# Run main function
main
