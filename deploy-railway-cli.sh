#!/bin/bash

# ==============================================================================
# CriOS Dr. Crowe Coder - Railway CLI Deployment Script
# ==============================================================================
#
# Run this script from your LOCAL MACHINE (not in the sandboxed environment)
#
# Prerequisites:
# - Node.js 18+ installed (check: node -v)
# - Git repository cloned
# - Railway account created (https://railway.app)
# - Anthropic API key ready (https://console.anthropic.com)
#
# Usage:
#   chmod +x deploy-railway-cli.sh
#   ./deploy-railway-cli.sh
#
# ==============================================================================

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# ASCII Art
clear
echo -e "${CYAN}"
cat << "EOF"
   ____      _  ___  ____
  / ___|_ __(_)/ _ \/ ___|
 | |   | '__| | | | \___ \
 | |___| |  | | |_| |___) |
  \____|_|  |_|\___/|____/

  Dr. Crowe Coder - Railway Deployment
  AI-Powered Drug Discovery Platform
EOF
echo -e "${NC}"

# ==============================================================================
# Helper Functions
# ==============================================================================

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[✓]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[⚠]${NC} $1"
}

log_error() {
    echo -e "${RED}[✗]${NC} $1"
}

log_step() {
    echo ""
    echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${MAGENTA}  STEP $1: $2${NC}"
    echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
}

prompt_user() {
    echo -e "${CYAN}$1${NC}"
    read -p "> " response
    echo "$response"
}

confirm() {
    echo -e "${YELLOW}$1 (y/n)${NC}"
    read -p "> " -n 1 -r
    echo
    [[ $REPLY =~ ^[Yy]$ ]]
}

# ==============================================================================
# Pre-flight Checks
# ==============================================================================

log_step "1" "Pre-flight Checks"

# Check Node.js
log_info "Checking Node.js installation..."
if ! command -v node &> /dev/null; then
    log_error "Node.js is not installed!"
    echo "Please install Node.js 18+ from: https://nodejs.org"
    exit 1
fi
NODE_VERSION=$(node -v)
log_success "Node.js installed: $NODE_VERSION"

# Check npm
log_info "Checking npm installation..."
if ! command -v npm &> /dev/null; then
    log_error "npm is not installed!"
    exit 1
fi
NPM_VERSION=$(npm -v)
log_success "npm installed: $NPM_VERSION"

# Check git
log_info "Checking Git installation..."
if ! command -v git &> /dev/null; then
    log_error "Git is not installed!"
    exit 1
fi
GIT_VERSION=$(git --version)
log_success "Git installed: $GIT_VERSION"

# Check if we're in the right directory
log_info "Checking current directory..."
if [ ! -f "railway.toml" ]; then
    log_error "railway.toml not found!"
    log_warning "Please run this script from the crios-dr-crowe-coder repository root"
    echo ""
    echo "Try:"
    echo "  cd /path/to/crios-dr-crowe-coder"
    echo "  git checkout claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH"
    echo "  ./deploy-railway-cli.sh"
    exit 1
fi
log_success "Found railway.toml - in correct directory"

# Check if on correct branch
CURRENT_BRANCH=$(git branch --show-current)
if [ "$CURRENT_BRANCH" != "claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH" ]; then
    log_warning "Current branch: $CURRENT_BRANCH"
    log_warning "Expected branch: claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH"

    if confirm "Switch to correct branch?"; then
        git checkout claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH
        git pull origin claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH
        log_success "Switched to correct branch"
    else
        log_warning "Continuing on current branch..."
    fi
fi

log_success "All pre-flight checks passed!"

# ==============================================================================
# Install Railway CLI
# ==============================================================================

log_step "2" "Install Railway CLI"

if command -v railway &> /dev/null; then
    RAILWAY_VERSION=$(railway --version 2>&1)
    log_success "Railway CLI already installed: $RAILWAY_VERSION"

    if confirm "Update to latest version?"; then
        log_info "Updating Railway CLI..."
        npm update -g @railway/cli
        log_success "Railway CLI updated"
    fi
else
    log_info "Installing Railway CLI..."
    echo "This may take 1-2 minutes..."
    npm i -g @railway/cli

    if command -v railway &> /dev/null; then
        log_success "Railway CLI installed successfully!"
    else
        log_error "Railway CLI installation failed!"
        echo "Please install manually:"
        echo "  npm i -g @railway/cli"
        exit 1
    fi
fi

# ==============================================================================
# Railway Authentication
# ==============================================================================

log_step "3" "Railway Authentication"

log_info "Checking Railway authentication..."

if railway whoami &> /dev/null; then
    RAILWAY_USER=$(railway whoami 2>&1)
    log_success "Already logged in as: $RAILWAY_USER"

    if ! confirm "Continue with this account?"; then
        log_info "Logging out..."
        # Railway CLI doesn't have logout, so we'll just re-login
        railway login
    fi
else
    log_info "Not logged in to Railway"
    log_warning "This will open your browser for authentication"

    if confirm "Ready to login to Railway?"; then
        railway login

        if railway whoami &> /dev/null; then
            RAILWAY_USER=$(railway whoami 2>&1)
            log_success "Logged in as: $RAILWAY_USER"
        else
            log_error "Railway login failed!"
            exit 1
        fi
    else
        log_error "Railway authentication required to continue"
        exit 1
    fi
fi

# ==============================================================================
# Project Setup
# ==============================================================================

log_step "4" "Railway Project Setup"

if [ -d ".railway" ]; then
    log_success "Railway project already linked"

    if confirm "Use existing project?"; then
        PROJECT_INFO=$(railway status 2>&1 | head -n 5 || echo "Project linked")
        echo "$PROJECT_INFO"
    else
        log_info "Creating new project..."
        rm -rf .railway
        railway init
    fi
else
    log_info "No Railway project found"

    if confirm "Create new Railway project?"; then
        railway init
        log_success "Railway project created"
    else
        log_warning "You can also link an existing project:"
        echo "  railway link"

        if confirm "Link existing project now?"; then
            railway link
        else
            log_error "Project setup required to continue"
            exit 1
        fi
    fi
fi

# ==============================================================================
# Database Provisioning
# ==============================================================================

log_step "5" "Database Provisioning"

# Check for PostgreSQL
log_info "Checking for PostgreSQL..."
if railway variables 2>&1 | grep -q "POSTGRES_HOST"; then
    log_success "PostgreSQL already provisioned"
else
    log_info "PostgreSQL not found"

    if confirm "Add PostgreSQL database?"; then
        log_info "Provisioning PostgreSQL... (this takes ~1 minute)"
        railway add --plugin postgresql

        # Wait for provisioning
        sleep 10

        if railway variables 2>&1 | grep -q "POSTGRES_HOST"; then
            log_success "PostgreSQL provisioned successfully!"
        else
            log_warning "PostgreSQL may still be provisioning..."
            log_warning "Check Railway dashboard if deployment fails"
        fi
    else
        log_warning "Skipping PostgreSQL - deployment may fail"
    fi
fi

# Check for Redis
log_info "Checking for Redis..."
if railway variables 2>&1 | grep -q "REDIS_HOST"; then
    log_success "Redis already provisioned"
else
    log_info "Redis not found"

    if confirm "Add Redis cache?"; then
        log_info "Provisioning Redis... (this takes ~1 minute)"
        railway add --plugin redis

        # Wait for provisioning
        sleep 10

        if railway variables 2>&1 | grep -q "REDIS_HOST"; then
            log_success "Redis provisioned successfully!"
        else
            log_warning "Redis may still be provisioning..."
            log_warning "Check Railway dashboard if deployment fails"
        fi
    else
        log_warning "Skipping Redis - worker may not function"
    fi
fi

# ==============================================================================
# Environment Variables
# ==============================================================================

log_step "6" "Environment Variables"

# Anthropic API Key
log_info "Checking ANTHROPIC_API_KEY..."
if railway variables 2>&1 | grep -q "ANTHROPIC_API_KEY"; then
    log_success "ANTHROPIC_API_KEY is already set"

    if confirm "Update ANTHROPIC_API_KEY?"; then
        API_KEY=$(prompt_user "Enter your Anthropic API key:")
        railway variables set ANTHROPIC_API_KEY="$API_KEY"
        log_success "ANTHROPIC_API_KEY updated"
    fi
else
    log_warning "ANTHROPIC_API_KEY not set (required!)"
    echo ""
    echo "Get your API key from: https://console.anthropic.com"
    echo ""

    API_KEY=$(prompt_user "Enter your Anthropic API key (or press Enter to skip):")

    if [ ! -z "$API_KEY" ]; then
        railway variables set ANTHROPIC_API_KEY="$API_KEY"
        log_success "ANTHROPIC_API_KEY set"
    else
        log_error "ANTHROPIC_API_KEY is required for deployment!"
        log_warning "You can set it later with:"
        echo "  railway variables set ANTHROPIC_API_KEY=your_key"

        if ! confirm "Continue without API key?"; then
            exit 1
        fi
    fi
fi

# JWT Secret
log_info "Checking JWT_SECRET..."
if railway variables 2>&1 | grep -q "JWT_SECRET"; then
    log_success "JWT_SECRET is already set"
else
    log_info "Generating JWT_SECRET..."
    JWT_SECRET=$(openssl rand -base64 32)
    railway variables set JWT_SECRET="$JWT_SECRET"
    log_success "JWT_SECRET generated and set"
fi

# Additional environment variables
log_info "Setting production environment variables..."
railway variables set ENVIRONMENT=production 2>&1 > /dev/null || true
railway variables set DEBUG=false 2>&1 > /dev/null || true
railway variables set LOG_LEVEL=info 2>&1 > /dev/null || true
railway variables set PYTHONUNBUFFERED=1 2>&1 > /dev/null || true
log_success "Production variables configured"

# ==============================================================================
# Deployment
# ==============================================================================

log_step "7" "Deploy to Railway"

log_info "Deployment configuration:"
echo "  - Backend API (FastAPI + RDKit)"
echo "  - Frontend (Next.js 14)"
echo "  - Celery Worker (ML Pipeline)"
echo "  - Core API (Discovery Engine)"
echo "  - PostgreSQL Database"
echo "  - Redis Cache"
echo ""
log_warning "First deployment takes ~10-15 minutes (installing RDKit, PyTorch, etc.)"
log_info "Subsequent deployments: ~3-5 minutes (cached layers)"
echo ""

if confirm "Start deployment now?"; then
    log_info "Deploying to Railway..."
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  Railway will now build and deploy all services"
    echo "  This will take approximately 10-15 minutes"
    echo "  You can monitor progress in the Railway dashboard"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""

    # Deploy with Railway
    railway up

    log_success "Deployment initiated!"
else
    log_warning "Deployment skipped"
    echo ""
    echo "Deploy manually with:"
    echo "  railway up"
    echo ""
    exit 0
fi

# ==============================================================================
# Post-Deployment
# ==============================================================================

log_step "8" "Post-Deployment Verification"

log_info "Waiting for services to start... (60 seconds)"
sleep 60

log_info "Checking service status..."
railway status

echo ""
log_info "Fetching service URLs..."
echo ""

# Try to get URLs
if command -v jq &> /dev/null; then
    # If jq is available, parse JSON
    BACKEND_URL=$(railway status --json 2>/dev/null | jq -r '.services[] | select(.name | contains("backend")) | .url' | head -n 1)
    FRONTEND_URL=$(railway status --json 2>/dev/null | jq -r '.services[] | select(.name | contains("frontend")) | .url' | head -n 1)

    if [ ! -z "$BACKEND_URL" ] && [ "$BACKEND_URL" != "null" ]; then
        echo -e "${GREEN}Backend URL:${NC} $BACKEND_URL"
    fi

    if [ ! -z "$FRONTEND_URL" ] && [ "$FRONTEND_URL" != "null" ]; then
        echo -e "${GREEN}Frontend URL:${NC} $FRONTEND_URL"
    fi
else
    log_info "Install jq to see service URLs automatically"
    log_info "Or view them in Railway dashboard"
fi

# ==============================================================================
# Health Checks
# ==============================================================================

log_step "9" "Health Checks"

if [ ! -z "$BACKEND_URL" ] && [ "$BACKEND_URL" != "null" ]; then
    log_info "Testing backend health endpoint..."

    if curl -f -s "${BACKEND_URL}/health" > /dev/null 2>&1; then
        log_success "Backend is healthy! ✓"
    else
        log_warning "Backend health check failed (may still be starting)"
        log_info "Check logs: railway logs --service crios-backend"
    fi
fi

if [ ! -z "$FRONTEND_URL" ] && [ "$FRONTEND_URL" != "null" ]; then
    log_info "Testing frontend accessibility..."

    if curl -f -s "$FRONTEND_URL" > /dev/null 2>&1; then
        log_success "Frontend is accessible! ✓"
    else
        log_warning "Frontend check failed (may still be building)"
        log_info "Check logs: railway logs --service crios-frontend"
    fi
fi

# ==============================================================================
# Deployment Summary
# ==============================================================================

echo ""
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${CYAN}              DEPLOYMENT COMPLETE! 🎉                   ${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
log_success "CriOS Dr. Crowe Coder deployed to Railway!"
echo ""

echo -e "${BLUE}Next Steps:${NC}"
echo ""
echo "  1. View your services:"
echo "     ${CYAN}railway status${NC}"
echo ""
echo "  2. Open Railway dashboard:"
echo "     ${CYAN}railway open${NC}"
echo ""
echo "  3. View backend logs:"
echo "     ${CYAN}railway logs --service crios-backend --follow${NC}"
echo ""
echo "  4. View frontend logs:"
echo "     ${CYAN}railway logs --service crios-frontend --follow${NC}"
echo ""
echo "  5. Test your deployment:"
if [ ! -z "$BACKEND_URL" ] && [ "$BACKEND_URL" != "null" ]; then
    echo "     ${CYAN}curl ${BACKEND_URL}/health${NC}"
fi
if [ ! -z "$FRONTEND_URL" ] && [ "$FRONTEND_URL" != "null" ]; then
    echo "     ${CYAN}open ${FRONTEND_URL}${NC}"
fi
echo ""

echo -e "${BLUE}Useful Commands:${NC}"
echo ""
echo "  View all environment variables:"
echo "    ${CYAN}railway variables${NC}"
echo ""
echo "  SSH into service:"
echo "    ${CYAN}railway run --service crios-backend bash${NC}"
echo ""
echo "  Restart service:"
echo "    ${CYAN}railway restart --service crios-backend${NC}"
echo ""
echo "  Add custom domain:"
echo "    ${CYAN}railway domain add yourdomain.com${NC}"
echo ""

echo -e "${BLUE}Documentation:${NC}"
echo "  • Complete Guide: RAILWAY_DEPLOYMENT.md"
echo "  • Quick Start: QUICK_START_RAILWAY.md"
echo "  • Troubleshooting: See RAILWAY_DEPLOYMENT.md"
echo ""

log_success "Happy discovering! 💊🔬"
echo ""
