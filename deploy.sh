#!/bin/bash
# CriOS Complete Deployment Script
# Deploys both frontend and backend to Vercel and Fly.io

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}"
echo "╔════════════════════════════════════════════════════════╗"
echo "║     CriOS Dr. Crowe Coder - Deployment Script         ║"
echo "║                                                        ║"
echo "║  Deploys:                                             ║"
echo "║    • Frontend to Vercel                               ║"
echo "║    • Algorithm Studio UI to Vercel                    ║"
echo "║    • Backend to Fly.io (optional)                     ║"
echo "╚════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# Check if Vercel CLI is installed
if ! command -v vercel &> /dev/null; then
    echo -e "${YELLOW}⚠️  Vercel CLI not found. Installing...${NC}"
    npm install -g vercel
fi

# Check if user is logged in to Vercel
echo -e "${BLUE}🔐 Checking Vercel authentication...${NC}"
if ! vercel whoami &> /dev/null; then
    echo -e "${YELLOW}⚠️  Not logged in to Vercel. Please login:${NC}"
    vercel login
    echo -e "${GREEN}✅ Logged in to Vercel${NC}"
else
    echo -e "${GREEN}✅ Already logged in to Vercel${NC}"
fi

# Ask what to deploy
echo ""
echo "What would you like to deploy?"
echo "1) Frontend only (platform/frontend)"
echo "2) Algorithm Studio UI only (ui)"
echo "3) Both Frontend and UI"
echo "4) Everything (Frontend, UI, and Backend)"
echo ""
read -p "Enter choice [1-4]: " deploy_choice

deploy_frontend() {
    echo ""
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}📦 Deploying Frontend to Vercel...${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

    cd platform/frontend

    # Check if npm dependencies are installed
    if [ ! -d "node_modules" ]; then
        echo -e "${YELLOW}📥 Installing dependencies...${NC}"
        npm install
    fi

    # Ask for environment variables
    echo ""
    echo -e "${YELLOW}Environment Variables:${NC}"
    read -p "Enter NEXT_PUBLIC_API_URL (default: https://crios-backend.fly.dev): " api_url
    api_url=${api_url:-https://crios-backend.fly.dev}

    read -p "Enter NEXT_PUBLIC_WS_URL (default: wss://crios-backend.fly.dev): " ws_url
    ws_url=${ws_url:-wss://crios-backend.fly.dev}

    # Create .env.production file
    cat > .env.production <<EOF
NEXT_PUBLIC_API_URL=$api_url
NEXT_PUBLIC_WS_URL=$ws_url
NODE_ENV=production
EOF

    echo -e "${GREEN}✅ Environment variables configured${NC}"

    # Deploy
    read -p "Deploy to production? (y/n): " prod
    if [ "$prod" = "y" ] || [ "$prod" = "Y" ]; then
        echo -e "${BLUE}🚀 Deploying to production...${NC}"
        vercel --prod --yes
    else
        echo -e "${BLUE}🚀 Deploying to preview...${NC}"
        vercel --yes
    fi

    cd ../..
    echo -e "${GREEN}✅ Frontend deployed successfully!${NC}"
}

deploy_ui() {
    echo ""
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}📦 Deploying Algorithm Studio UI to Vercel...${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

    cd ui

    # Check if npm dependencies are installed
    if [ ! -d "node_modules" ]; then
        echo -e "${YELLOW}📥 Installing dependencies...${NC}"
        npm install
    fi

    # Ask for environment variables
    echo ""
    echo -e "${YELLOW}Environment Variables:${NC}"
    read -p "Enter NEXT_PUBLIC_API_URL (default: https://crios-backend.fly.dev): " api_url
    api_url=${api_url:-https://crios-backend.fly.dev}

    # Create .env.production file
    cat > .env.production <<EOF
NEXT_PUBLIC_API_URL=$api_url
NODE_ENV=production
EOF

    echo -e "${GREEN}✅ Environment variables configured${NC}"

    # Deploy
    read -p "Deploy to production? (y/n): " prod
    if [ "$prod" = "y" ] || [ "$prod" = "Y" ]; then
        echo -e "${BLUE}🚀 Deploying to production...${NC}"
        vercel --prod --yes
    else
        echo -e "${BLUE}🚀 Deploying to preview...${NC}"
        vercel --yes
    fi

    cd ..
    echo -e "${GREEN}✅ Algorithm Studio UI deployed successfully!${NC}"
}

deploy_backend() {
    echo ""
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}📦 Deploying Backend to Fly.io...${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

    # Check if flyctl is installed
    if ! command -v flyctl &> /dev/null; then
        echo -e "${YELLOW}⚠️  Fly.io CLI not found. Installing...${NC}"
        curl -L https://fly.io/install.sh | sh
        export PATH="$HOME/.fly/bin:$PATH"
    fi

    # Check if user is logged in to Fly.io
    echo -e "${BLUE}🔐 Checking Fly.io authentication...${NC}"
    if ! flyctl auth whoami &> /dev/null; then
        echo -e "${YELLOW}⚠️  Not logged in to Fly.io. Please login:${NC}"
        flyctl auth login
        echo -e "${GREEN}✅ Logged in to Fly.io${NC}"
    else
        echo -e "${GREEN}✅ Already logged in to Fly.io${NC}"
    fi

    cd platform/backend

    # Check if app exists
    echo -e "${BLUE}📋 Checking if app exists...${NC}"
    if ! flyctl status -a crios-backend &> /dev/null; then
        echo -e "${YELLOW}Creating new Fly.io app...${NC}"
        flyctl apps create crios-backend

        # Ask for secrets
        echo ""
        echo -e "${YELLOW}Please provide the following secrets:${NC}"
        read -p "DATABASE_URL: " db_url
        read -p "REDIS_URL: " redis_url
        read -p "ANTHROPIC_API_KEY: " api_key

        flyctl secrets set \
            DATABASE_URL="$db_url" \
            REDIS_URL="$redis_url" \
            ANTHROPIC_API_KEY="$api_key" \
            -a crios-backend
    fi

    # Deploy
    echo -e "${BLUE}🚀 Deploying backend...${NC}"
    flyctl deploy

    cd ../..
    echo -e "${GREEN}✅ Backend deployed successfully!${NC}"
}

# Execute based on choice
case $deploy_choice in
    1)
        deploy_frontend
        ;;
    2)
        deploy_ui
        ;;
    3)
        deploy_frontend
        deploy_ui
        ;;
    4)
        deploy_frontend
        deploy_ui
        deploy_backend
        ;;
    *)
        echo -e "${RED}❌ Invalid choice. Exiting.${NC}"
        exit 1
        ;;
esac

# Summary
echo ""
echo -e "${GREEN}"
echo "╔════════════════════════════════════════════════════════╗"
echo "║              🎉 Deployment Complete! 🎉                ║"
echo "╚════════════════════════════════════════════════════════╝"
echo -e "${NC}"

echo -e "${BLUE}Your applications are now live!${NC}"
echo ""

if [ "$deploy_choice" = "1" ] || [ "$deploy_choice" = "3" ] || [ "$deploy_choice" = "4" ]; then
    echo -e "📱 Frontend: Check Vercel dashboard for URL"
fi

if [ "$deploy_choice" = "2" ] || [ "$deploy_choice" = "3" ] || [ "$deploy_choice" = "4" ]; then
    echo -e "🎨 Algorithm Studio: Check Vercel dashboard for URL"
fi

if [ "$deploy_choice" = "4" ]; then
    echo -e "⚙️  Backend: https://crios-backend.fly.dev"
    echo -e "📚 API Docs: https://crios-backend.fly.dev/docs"
fi

echo ""
echo -e "${YELLOW}Next steps:${NC}"
echo "1. Visit Vercel dashboard: https://vercel.com/dashboard"
echo "2. Test your deployments"
echo "3. Configure custom domains (optional)"
echo "4. Set up monitoring and alerts"
echo ""
echo -e "${GREEN}Documentation:${NC}"
echo "  • Quick Guide: DEPLOY_NOW.md"
echo "  • Full Guide: DEPLOYMENT.md"
echo ""
